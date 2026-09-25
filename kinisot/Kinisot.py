#!/usr/bin/python
"""Kinisot: kinetic and equilibrium isotope effects from computed Hessians.

Comments and/or additions are welcome (send e-mail to
robert.paton@colostate.edu).
"""

import os
import sys
import time
import warnings
from argparse import SUPPRESS, ArgumentParser

import numpy as np

from . import __version__
from .exceptions import KinisotError, KinisotInputError, KinisotWarning
from .Hess_to_Freq import is_linear, level_of_theory, mass_weight, parse_gaussian, read_hess, substitute
from .vib_scale_factors import scaling_data, scaling_refs

__all__ = [
    "compute_isotope_effect",
    "calc_rpfr",
    "harmonic_frequencies",
    "find_scaling_factor",
    "get_frequency_scaling",
    "read_hess",
    "is_linear",
    "level_of_theory",
    "Logger",
    "build_parser",
    "main",
    "__version__",
]

# PHYSICAL CONSTANTS (CODATA 2010; SI apart from the speed of light in cm/s)
PLANCK_CONSTANT = 6.62606957e-34  # J s
BOLTZMANN_CONSTANT = 1.3806488e-23  # J / K
SPEED_OF_LIGHT = 2.99792458e10  # cm / s
ENERGY_AU = 4.35974434e-18  # J
BOHR_RADIUS = 5.2917721092e-11  # m
ATOMIC_MASS_UNIT = 1.660538921e-27  # kg

# Multiply a mass-weighted Hessian in Hartree/(amu Bohr^2) by this to get eigenvalues in cm^-2
HESSIAN_TO_WAVENUMBER_SQ = ENERGY_AU / (BOHR_RADIUS**2 * ATOMIC_MASS_UNIT) / ((SPEED_OF_LIGHT * 2 * np.pi) ** 2)

# print formatting
space = "   "
dash = "--"
dash_line = space * 17 + " " + dash * 37


class Logger:
    """Writes results to the terminal and to a results file.

    New results are appended to ``path`` so that a series of runs (one per
    substitution, as in the bundled examples) accumulates in one file; pass
    ``overwrite=True`` to start afresh. Use as a context manager so the file
    is always closed.
    """

    def __init__(self, path="Kinisot_output.dat", quiet=False, overwrite=False):
        self.path = path
        self.quiet = quiet
        self.log = open(path, "w" if overwrite else "a")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.Finalize()
        return False

    def Write(self, message):
        """Write a message to the terminal (unless quiet) and to the file."""
        if not self.quiet:
            print(message, end="")
        self.log.write(message)

    def Writeonlyfile(self, message):
        """Write a message only to the file and not to the terminal."""
        self.log.write("\n" + message + "\n")

    def Finalize(self):
        if not self.log.closed:
            self.log.close()


def find_scaling_factor(level):
    """
    Look up the ZPE vibrational scaling factor for a level of theory.

    The level string (as written in the Gaussian archive, e.g. RM062X/MG3S)
    is matched exactly against the Truhlar database entries after normalizing
    case and hyphens and stripping Gaussian's R/U/RO spin prefix. Returns
    (factor, reference) or (None, None) if the level is not in the database.
    """

    def norm(name):
        return name.upper().replace("-", "")

    candidates = {norm(level)}
    for prefix in ("RO", "R", "U"):
        if level.upper().startswith(prefix):
            candidates.add(norm(level[len(prefix) :]))
    for scal in scaling_data:
        if norm(scal["level"].decode("utf-8")) in candidates:
            return scal["zpe_fac"], scaling_refs[scal["zpe_ref"]]
    return None, None


def get_frequency_scaling(files, log):
    """Detect the level of theory of all files and look up the ZPE scaling factor.

    Writes what was found to the log and returns the factor (1.0 when the
    files disagree or the level of theory is not in the database).
    """
    levels = {}
    for file in files:
        levels.setdefault(level_of_theory(file) or "unknown", []).append(file)

    if len(levels) > 1:
        log.Write("\n  WARNING: the files were not computed at the same level of theory:")
        for level, names in levels.items():
            log.Write("\n     %s: %s" % (level, ", ".join(names)))
        log.Write("\n  Unable to assign a vibrational scaling factor; using 1.0 (override with -s)")
        return 1.0

    level = next(iter(levels))
    factor, ref = find_scaling_factor(level) if level != "unknown" else (None, None)
    if factor is None:
        log.Write("\n  Unable to find vibrational scaling factor for %s; using value of 1.0 (override with -s)" % level)
        return 1.0
    factor = round(float(factor), 4)
    log.Write(
        "\n  Found vibrational scaling factor %s for %s level of theory "
        "(ZPE factor, Truhlar database)" % (factor, level)
    )
    log.Write("\n  REF: " + ref)
    return factor


def _reduced_energies(frequency_wn, temperature):
    """u = h c nu / k T for an array of wavenumbers."""
    return PLANCK_CONSTANT * SPEED_OF_LIGHT * np.asarray(frequency_wn, dtype=float) / (BOLTZMANN_CONSTANT * temperature)


def calc_product_factor(frequency_wn):
    """
    Log of the product of (scaled) vibrational frequencies, expressed as
    vibrational temperatures, which gives the Teller-Redlich product factor.
    There is no temperature dependence to this term in the BM equation.
    Everything is done logarithmically to avoid big numbers.
    """
    hv_over_k = PLANCK_CONSTANT * SPEED_OF_LIGHT * np.asarray(frequency_wn, dtype=float) / BOLTZMANN_CONSTANT
    return float(np.sum(np.log(hv_over_k)))


def calc_zpe_factor(frequency_wn, temperature):
    """
    Log of the ZPE term of the BM equation, sum(u/2), for (scaled) vibrational
    frequencies. ZPEs themselves are not temperature dependent although the
    exponential form of this term is.
    """
    return float(0.5 * np.sum(_reduced_energies(frequency_wn, temperature)))


def calc_excitation_factor(frequency_wn, temperature):
    """
    Log of the excitation term of the BM equation, sum(ln(1 - exp(-u))), for
    (scaled) vibrational frequencies. This term is temperature dependent.
    """
    u = _reduced_energies(frequency_wn, temperature)
    return float(np.sum(np.log1p(-np.exp(-u))))


def harmonic_frequencies(mw_hessian):
    """Harmonic frequencies in cm-1 (negative for imaginary modes), ascending."""
    eigenvalues = np.linalg.eigvalsh(mw_hessian * HESSIAN_TO_WAVENUMBER_SQ)
    return np.copysign(np.sqrt(np.abs(eigenvalues)), eigenvalues)


class calc_rpfr:
    """Reduced isotopic partition function ratio terms for one side of a reaction.

    ``files`` are one or more frequency outputs (several for bimolecular
    reactions) and ``isomer`` the matching isotope labels. The mass-weighted
    Hessian of each file is diagonalized, the 5/6 external modes (and one
    reaction-coordinate mode, if present) are removed, and the logarithmic
    ZPE, excitation and product terms are accumulated in ``ZPE``, ``EXC``
    and ``PF``.

    Attributes set per instance: ``frequency_wn`` (kept modes, all files),
    ``discarded_wn`` (external modes removed, per file), ``im_frequencies``
    (magnitude of the reaction-coordinate mode, per file that has one),
    ``im_frequency_wn`` (the same when exactly one file has one) and
    ``substituted`` (Substitution records).
    """

    def __init__(self, files, isomer, temperature=298.15, freq_scale_factor=1.0, freq_cutoff=50.0):
        self.files = list(files)
        self.isomer = list(isomer)
        self.PF, self.ZPE, self.EXC = 0.0, 0.0, 0.0
        self.frequency_wn = []
        self.kept_wn = {}
        self.discarded_wn = {}
        self.im_frequencies = {}
        self.substituted = []

        for file, iso in zip(self.files, self.isomer):
            data = parse_gaussian(file)
            masses, applied = substitute(data, iso)
            self.substituted.extend(applied)

            # Frequencies (scaled) from the mass-weighted Hessian in Hartree/(amu Bohr^2)
            freqs = harmonic_frequencies(mass_weight(data.hessian, masses)) * freq_scale_factor

            # 5 or 6 external modes are removed (linear / non-linear molecule),
            # plus one reaction-coordinate mode when it is imaginary beyond the cutoff
            n_external = 5 if data.linear else 6
            imaginary = freqs[freqs < -freq_cutoff]
            if len(imaginary):
                self.im_frequencies[file] = -float(freqs[0])
                n_external += 1
                if len(imaginary) > 1:
                    warnings.warn(
                        "%s has %d imaginary frequencies beyond the %.1f cm-1 cutoff (%s); only "
                        "the largest is treated as the reaction coordinate and the others are "
                        "discarded with the external modes, which leaves one low-frequency "
                        "external mode in the vibrational product. Check the structure."
                        % (file, len(imaginary), freq_cutoff, ", ".join("%.1fi" % -f for f in imaginary)),
                        KinisotWarning,
                        stacklevel=2,
                    )
            discarded = freqs[(1 if len(imaginary) else 0) : n_external]
            kept = freqs[n_external:]
            if np.any(kept <= 0):
                raise KinisotInputError(
                    "%s: %d non-positive frequencies remain after removing %d external modes "
                    "(%s); the structure is not a stationary point Kinisot can use"
                    % (file, int(np.sum(kept <= 0)), n_external, ", ".join("%.1f" % f for f in kept[kept <= 0]))
                )
            self.discarded_wn[file] = [float(f) for f in discarded]
            self.kept_wn[file] = [float(f) for f in kept]
            self.frequency_wn.extend(self.kept_wn[file])

            # Calculate the excitation factor (EXC), the ZPE (ZPE) and Teller-Redlich product factor (PF)
            self.PF += calc_product_factor(kept)
            self.ZPE += calc_zpe_factor(kept, temperature)
            self.EXC += calc_excitation_factor(kept, temperature)

        if len(self.im_frequencies) > 1:
            raise KinisotInputError(
                "more than one file on the same side of the reaction has an imaginary "
                "frequency: %s" % ", ".join("%s (%.1fi cm-1)" % item for item in self.im_frequencies.items())
            )
        if len(self.im_frequencies) == 1:
            self.im_frequency_wn = next(iter(self.im_frequencies.values()))


def _describe(substitutions):
    return ", ".join(str(s) for s in substitutions) if substitutions else "none"


def _check_substitution_balance(reactant_side, other_side, side_name):
    """Both sides must carry the same isotopic substitutions (same elements)."""
    left = sorted(s.symbol for s in reactant_side.substituted)
    right = sorted(s.symbol for s in other_side.substituted)
    if not left and not right:
        raise KinisotInputError("no isotopic substitution requested: every --iso label is '0'")
    if left != right:
        raise KinisotInputError(
            "the isotopic substitutions differ between the reactant side [%s] and the %s side "
            "[%s]. An isotope effect compares the same isotopologue on both sides; check the atom "
            "numbering in each file."
            % (_describe(reactant_side.substituted), side_name, _describe(other_side.substituted))
        )


def compute_isotope_effect(rct, ts, prd, label, temperature=298.15, freq_scale_factor=1.0, freq_cutoff=50.0):
    """Compute a KIE (``ts`` given) or EQE (``prd`` given) from Gaussian outputs.

    ``rct`` and ``ts``/``prd`` are lists of files; ``label`` holds one isotope
    label per file, reactant files first ('0' for a file without substitution).

    Returns ``(species, ZPE, EXC, TRPF, KIE, KIE_tunnel, tunnel_corr, freq_fac)``
    where ``species`` lists the four calc_rpfr results (reactant light/heavy,
    TS-or-product light/heavy). Raises KinisotInputError for inconsistent
    input and KinisotParseError for unusable files.
    """
    if (ts is None) == (prd is None):
        raise KinisotInputError(
            "give either transition structure files (KIE) or product files (EQE), not both and not neither"
        )
    rct = list(rct)
    other = list(ts if ts is not None else prd)
    side_name = "transition structure" if ts is not None else "product"
    label = list(label)
    if len(label) != len(rct) + len(other):
        raise KinisotInputError(
            "%d files were given but %d isotope labels; give one label per file" % (len(rct) + len(other), len(label))
        )
    if temperature <= 0:
        raise KinisotInputError("temperature must be positive (got %s K)" % temperature)
    if freq_scale_factor <= 0:
        raise KinisotInputError("the vibrational scaling factor must be positive (got %s)" % freq_scale_factor)

    # Calculates the RPFR terms for each species and its isotopologue
    species = []
    for iso in [["0"] * len(rct), label[0 : len(rct)]]:
        species.append(calc_rpfr(rct, iso, temperature, freq_scale_factor, freq_cutoff))
    for iso in [["0"] * len(other), label[len(rct) :]]:
        species.append(calc_rpfr(other, iso, temperature, freq_scale_factor, freq_cutoff))

    # Reactants (and products) must be minima
    minima = [("reactant", species[0])] + ([("product", species[2])] if prd is not None else [])
    for role, sp in minima:
        if sp.im_frequencies:
            file, im = next(iter(sp.im_frequencies.items()))
            raise KinisotInputError(
                "%s has an imaginary frequency (%.1fi cm-1) but was given as a %s. Reactants and "
                "products must be minima; re-optimize the structure or, if this is a spurious "
                "low mode, raise --imag-cutoff above %.1f" % (file, im, role, im)
            )

    if ts is not None:
        # Check for the presence of an imaginary frequency in both TS isotopologues
        if not (species[2].im_frequencies and species[3].im_frequencies):
            raise KinisotInputError(
                "Kinisot requires a transition structure with an imaginary frequency beyond the "
                "%.1f cm-1 cutoff (--imag-cutoff)! Files given as --ts: %s" % (freq_cutoff, ", ".join(other))
            )
        freq_fac = species[2].im_frequency_wn / species[3].im_frequency_wn
    else:
        freq_fac = 1.0

    _check_substitution_balance(species[1], species[3], side_name)

    # Application of the Bigeleisen-Mayer equation
    ZPE = np.e ** (species[0].ZPE - species[1].ZPE - species[2].ZPE + species[3].ZPE)
    EXC = np.e ** (species[0].EXC - species[1].EXC - species[2].EXC + species[3].EXC)
    TRPF = np.e ** (species[2].PF - species[3].PF - species[0].PF + species[1].PF)

    # A correction factor for QM-tunneling (Bell infinite parabola)
    # Conversion from wavenumbers to SI energy units; then divide by kT
    tofreq = SPEED_OF_LIGHT * PLANCK_CONSTANT / BOLTZMANN_CONSTANT / temperature
    if ts is not None:
        parabolic_tunn_corr = (
            freq_fac
            * np.sin(0.5 * tofreq * species[3].im_frequency_wn)
            / np.sin(0.5 * tofreq * species[2].im_frequency_wn)
        )
    else:
        parabolic_tunn_corr = 1.0

    # (a) the Bigeleisen-Mayer KIE with classical nuclei and
    # (b) a value corrected to include quantum tunneling effects
    KIE_no_tunnel = freq_fac * ZPE * EXC * TRPF
    KIE_tunnel = KIE_no_tunnel * parabolic_tunn_corr

    return species, ZPE, EXC, TRPF, KIE_no_tunnel, KIE_tunnel, parabolic_tunn_corr, freq_fac


def _short_name(file):
    """File name without directory and extension, for the results table."""
    return os.path.splitext(os.path.basename(file))[0]


def _side_name(files):
    return " + ".join(_short_name(f) for f in files)


def write_results(log, files, labels, n_rct, is_kie, temperature, scale, result):
    """Write the results table for one calculation to the log."""
    species, ZPE, EXC, TRPF, KIE_no_tunnel, KIE_tunnel, parabolic_tunn_corr, freq_fac = result
    rct_files, oth_files = files[:n_rct], files[n_rct:]
    rct_name, oth_name = _side_name(rct_files), _side_name(oth_files)
    rct_iso, oth_iso = " / ".join(labels[:n_rct]), " / ".join(labels[n_rct:])

    log.Write("\n\n" + (space * 17) + "  Temp = " + str(temperature) + "K / Vib. scale factor = " + str(scale))
    log.Write(("\n  ").ljust(50))
    log.Write(
        " {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} \n".format(
            "V-ratio", "ZPE", "EXC", "TRPF", "KIE", "1D-tunn", "corr-KIE"
        )
    )

    # Per-species Bigeleisen-Mayer factors (light / heavy); the final line is their ratio
    def factors(light, heavy):
        return (np.e ** (light.ZPE - heavy.ZPE), np.e ** (light.EXC - heavy.EXC), np.e ** (heavy.PF - light.PF))

    log.Write("\no " + rct_name.ljust(47) + "   " + dash * 37)
    log.Write("\no " + oth_name.ljust(47))
    if is_kie:
        log.Write("{:10.1f}".format(species[2].im_frequency_wn))
    log.Write("\no " + (rct_name + ": iso @ " + rct_iso).ljust(47))
    log.Write("           {:10.3e} {:10.3e} {:10.3e}".format(*factors(species[0], species[1])))
    log.Write("\no " + (oth_name + ": iso @ " + oth_iso).ljust(47))
    if is_kie:
        log.Write(
            "{:10.1f} {:10.3e} {:10.3e} {:10.3e}".format(species[3].im_frequency_wn, *factors(species[2], species[3]))
        )
    else:
        log.Write("{:21.3e} {:10.3e} {:10.3e}".format(*factors(species[2], species[3])))

    log.Write("\n" + dash_line)
    log.Write(("\n  " + ("KIE" if is_kie else "EQE") + " @ " + str(temperature) + " K").ljust(50))
    if is_kie:
        log.Write(
            "{:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                freq_fac, ZPE, EXC, TRPF, KIE_no_tunnel, parabolic_tunn_corr, KIE_tunnel
            )
        )
    else:
        log.Write(
            "{:21.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                ZPE, EXC, TRPF, KIE_no_tunnel, parabolic_tunn_corr, KIE_tunnel
            )
        )
    log.Write("\n" + dash_line + "\n")

    # Which modes went into the partition functions, so that a misassigned external mode is visible
    log.Write("\n  Vibrational modes (scaled, cm-1): kept in the partition function / discarded as external modes\n")
    for sp, iso_labels in (
        (species[0], None),
        (species[1], labels[:n_rct]),
        (species[2], None),
        (species[3], labels[n_rct:]),
    ):
        for i, file in enumerate(sp.files):
            tag = "light" if iso_labels is None else "iso @ " + iso_labels[i]
            line = "  %s (%s):" % (_short_name(file), tag)
            if file in sp.im_frequencies:
                line += " imaginary %.1fi;" % sp.im_frequencies[file]
            line += " %d kept; discarded: %s" % (
                len(sp.kept_wn[file]),
                " ".join("%.1f" % f for f in sp.discarded_wn[file]),
            )
            log.Write(line + "\n")


def build_parser():
    parser = ArgumentParser(
        prog="kinisot",
        description="Kinetic (--ts) and equilibrium (--prd) isotope effects from Gaussian frequency "
        "calculations, using the Bigeleisen-Mayer equation and a Bell tunnelling correction.",
        epilog="Example: kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 5 -t 393 -s 0.961",
    )
    parser.add_argument(
        "--rct",
        dest="rct",
        action="append",
        required=True,
        metavar="FILE",
        help="reactant frequency output; repeat for bimolecular reactions",
    )
    parser.add_argument(
        "--ts", dest="ts", action="append", metavar="FILE", help="transition structure frequency output (KIE)"
    )
    parser.add_argument("--prd", dest="prd", action="append", metavar="FILE", help="product frequency output (EQE)")
    parser.add_argument(
        "--iso",
        dest="label",
        action="append",
        required=True,
        metavar="ATOMS",
        help="atom number(s) to replace with the heavy isotope (2H, 13C, 17O), comma "
        "separated, e.g. 7,8. Give one --iso per file in the order of the --rct "
        "then --ts/--prd files, or a single --iso when the atom numbering is the "
        "same in all files. Use 0 for a file without substitution.",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        dest="temperature",
        type=float,
        default=298.15,
        help="temperature in Kelvin (default 298.15)",
    )
    parser.add_argument(
        "-s",
        "--scale",
        dest="freq_scale_factor",
        type=float,
        default=None,
        help="vibrational scaling factor (default: ZPE factor from the Truhlar database "
        "for the detected level of theory, else 1.0)",
    )
    parser.add_argument(
        "--imag-cutoff",
        dest="freq_cutoff",
        type=float,
        default=50.0,
        help="a mode below -CUTOFF cm-1 is the reaction coordinate (default 50)",
    )
    parser.add_argument("--cutoff", dest="freq_cutoff", type=float, default=50.0, help=SUPPRESS)
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="Kinisot_output.dat",
        metavar="FILE",
        help="results file; new results are appended (default Kinisot_output.dat)",
    )
    parser.add_argument("--overwrite", action="store_true", help="start a fresh results file instead of appending")
    parser.add_argument("-q", "--quiet", action="store_true", help="do not print results to the terminal")
    parser.add_argument("--version", action="version", version="Kinisot " + __version__)
    return parser


def main(argv=None):
    """Command-line entry point. Returns the process exit code."""
    parser = build_parser()
    options = parser.parse_args(argv)

    if options.ts is None and options.prd is None:
        parser.error("either --ts (for a KIE) or --prd (for an EQE) is required")
    if options.ts is not None and options.prd is not None:
        parser.error("--ts and --prd cannot be combined: use --ts for a KIE or --prd for an EQE")
    if options.temperature <= 0:
        parser.error("the temperature must be positive")
    if options.freq_scale_factor is not None and options.freq_scale_factor <= 0:
        parser.error("the scaling factor must be positive")

    is_kie = options.ts is not None
    files = options.rct + (options.ts if is_kie else options.prd)
    # if only one set of labels is provided, assume that the atom numbering is the same for rct and ts/prd
    labels = options.label * 2 if len(options.label) == 1 else list(options.label)
    if len(labels) != len(files):
        parser.error(
            "%d files were given but %d --iso labels: give one --iso per file (0 for no "
            "substitution) or a single --iso when the atom numbering is the same" % (len(files), len(options.label))
        )

    try:
        log = Logger(options.output, quiet=options.quiet, overwrite=options.overwrite)
    except OSError as err:
        print(
            "\no  ERROR: cannot open the results file %s: %s\n" % (options.output, err.strerror or err), file=sys.stderr
        )
        return 1
    with log:
        log.Write("\n  " + "KINISOT.py v " + __version__ + ": " + time.strftime("%Y-%m-%d %H:%M") + "\n")
        for file, label in zip(files, labels):
            log.Write("  Species: {} isotopologue: {}\n".format(file, label))
        try:
            # if not specified try to automatically determine the vibrational scaling factor
            scale = options.freq_scale_factor
            if scale is None:
                scale = get_frequency_scaling(files, log)

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", KinisotWarning)
                result = compute_isotope_effect(
                    options.rct, options.ts, options.prd, labels, options.temperature, scale, options.freq_cutoff
                )
            for warning in caught:
                log.Write("\n  WARNING: " + str(warning.message))

            write_results(log, files, labels, len(options.rct), is_kie, options.temperature, scale, result)
        except KinisotError as err:
            log.Writeonlyfile("o  ERROR: " + str(err))
            print("\no  ERROR: " + str(err) + "\n", file=sys.stderr)
            return 1
        if not options.quiet:
            print("\n  Results appended to " + options.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
