"""Command-line interface: argument parsing, the results table and file outputs."""

import csv
import json
import os
import sys
import time
import warnings
from argparse import SUPPRESS, ArgumentParser

from . import __version__
from .api import compute_kie
from .ensemble import WEIGHTS, Conformers, EnsembleIsotopeEffect
from .exceptions import KinisotError, KinisotInputError, KinisotWarning
from .pathways import ChannelIsotopeEffect, SeriesIsotopeEffect

__all__ = [
    "Logger",
    "build_parser",
    "write_results",
    "write_ensemble_results",
    "write_series_results",
    "write_channel_results",
    "read_energies",
    "write_json",
    "append_csv",
    "main",
]

# print formatting
SPACE = "   "
DASH = "--"
DASH_LINE = SPACE * 17 + " " + DASH * 37


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


def parse_temperatures(text):
    """'393' -> [393.0]; '273,298,323' -> list; '250:350:10' -> range (inclusive end)."""
    values = []
    for part in str(text).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            pieces = part.split(":")
            if len(pieces) != 3:
                raise ValueError("temperature range must be START:STOP:STEP (got %r)" % part)
            start, stop, step = (float(p) for p in pieces)
            if step <= 0 or stop < start:
                raise ValueError("temperature range %r must have STOP >= START and STEP > 0" % part)
            t = start
            while t <= stop + 1e-9:
                values.append(round(t, 6))
                t += step
        else:
            values.append(float(part))
    if not values:
        raise ValueError("no temperature given")
    return values


def write_results(log, results):
    """Write the results table for one calculation (one or more temperatures) to the log."""
    results = list(results)
    result = results[0]
    is_kie = result.kind == "KIE"
    rct, oth = result.reactant, result.other
    rct_iso, oth_iso = " / ".join(rct.heavy.labels), " / ".join(oth.heavy.labels)

    log.Write(
        "\n\n"
        + (SPACE * 17)
        + "  Temp = "
        + (str(result.temperature) if len(results) == 1 else ", ".join(str(r.temperature) for r in results))
        + "K / Vib. scale factor = "
        + str(result.scale_factor)
    )
    if is_kie and result.tunneling != "bell":
        log.Write(" / tunnelling: " + result.tunneling)
        if result.barrier is not None:
            log.Write(" (barrier %.2f kcal/mol)" % result.barrier)
    if result.project:
        log.Write(" / external modes projected")
    log.Write(("\n  ").ljust(50))
    log.Write(
        " {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} \n".format(
            "V-ratio", "ZPE", "EXC", "TRPF", "KIE", "1D-tunn", "corr-KIE"
        )
    )

    # Per-species Bigeleisen-Mayer factors (light / heavy); the final line is their ratio
    log.Write("\no " + rct.name.ljust(47) + "   " + DASH * 37)
    log.Write("\no " + oth.name.ljust(47))
    if is_kie:
        log.Write("{:10.1f}".format(oth.light.imaginary))
    log.Write("\no " + (rct.name + ": iso @ " + rct_iso).ljust(47))
    log.Write("           {:10.3e} {:10.3e} {:10.3e}".format(rct.zpe_factor, rct.exc_factor, rct.trpf_factor))
    log.Write("\no " + (oth.name + ": iso @ " + oth_iso).ljust(47))
    if is_kie:
        log.Write(
            "{:10.1f} {:10.3e} {:10.3e} {:10.3e}".format(
                oth.heavy.imaginary, oth.zpe_factor, oth.exc_factor, oth.trpf_factor
            )
        )
    else:
        log.Write("{:21.3e} {:10.3e} {:10.3e}".format(oth.zpe_factor, oth.exc_factor, oth.trpf_factor))

    log.Write("\n" + DASH_LINE)
    for r in results:
        log.Write(("\n  " + r.kind + " @ " + str(r.temperature) + " K").ljust(50))
        if is_kie:
            log.Write(
                "{:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                    r.imag_ratio, r.zpe, r.exc, r.trpf, r.kie, r.tunnel_corr, r.kie_tunnel
                )
            )
        else:
            log.Write(
                "{:21.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                    r.zpe, r.exc, r.trpf, r.kie, r.tunnel_corr, r.kie_tunnel
                )
            )
        if r.reference is not None:
            ref = r.reference
            log.Write(
                ("\n  relative to iso @ " + " / ".join(ref.reactant.heavy.labels + ref.other.heavy.labels) + ":").ljust(
                    50
                )
            )
            log.Write(
                "{:>43} {:10.6f} {:>10} {:10.6f}".format(
                    "%s = %.6f" % (r.kind, ref.kie), r.kie_relative, "", r.kie_tunnel_relative
                )
                if is_kie
                else "{:>54} {:10.6f} {:>10} {:10.6f}".format(
                    "%s = %.6f" % (r.kind, ref.kie), r.kie_relative, "", r.kie_tunnel_relative
                )
            )
    log.Write("\n" + DASH_LINE + "\n")

    # Which modes went into the partition functions, so that a misassigned external mode is visible
    log.Write(
        "\n  Vibrational modes (scaled, cm-1): kept in the partition function / %s as external modes\n"
        % ("projected out (residual values)" if result.project else "discarded")
    )
    for isotopologue, tags in (
        (rct.light, None),
        (rct.heavy, rct.heavy.labels),
        (oth.light, None),
        (oth.heavy, oth.heavy.labels),
    ):
        for i, species in enumerate(isotopologue.species):
            tag = "light" if tags is None else "iso @ " + tags[i]
            line = "  %s (%s):" % (species.name, tag)
            if species.imaginary is not None:
                line += " imaginary %.1fi;" % species.imaginary
            line += " %d kept; discarded: %s" % (
                len(species.frequencies),
                " ".join("%.1f" % f for f in species.discarded),
            )
            log.Write(line + "\n")


def write_ensemble_results(log, results):
    """Write the conformer table and the ensemble isotope effects (one line per temperature) to the log."""
    results = list(results)
    result = results[0]
    is_kie = result.kind == "KIE"
    log.Write(
        "\n\n"
        + (SPACE * 17)
        + "  Temp = "
        + (str(result.temperature) if len(results) == 1 else ", ".join(str(r.temperature) for r in results))
        + "K / Vib. scale factor = "
        + str(result.scale_factor)
        + " / weights: "
        + result.weights
    )
    if is_kie and result.tunneling != "bell":
        log.Write(" / tunnelling: " + result.tunneling)
        if result.barrier is not None:
            log.Write(" (barrier %.2f kcal/mol)" % result.barrier)
    if result.project:
        log.Write(" / external modes projected")

    # the conformer table: each conformer's weight and its KIE against the ensemble on the other side
    log.Write(
        "\n\n  Conformers at %s K%s:"
        % (result.temperature, " (the JSON output has the other temperatures)" if len(results) > 1 else "")
    )
    log.Write(
        "\n  " + " ".join(["{:<44}".format(""), "{:>8} {:>5} {:>7} {:>10} {:>10} {:>10}"])
        .format("dG", "g", "pop %", "V-ratio", result.kind, "corr-" + result.kind if is_kie else "")
        .rstrip()
    )  # fmt: skip
    role = None
    for c in result.conformers:
        if (c.role, c.species) != role:
            role = (c.role, c.species)
            log.Write("\no %s %d, iso @ %s" % (c.role, c.species + 1, c.label))
        energy = "{:8.2f}".format(c.free_energy) if c.free_energy is not None else "{:>8}".format("-")
        ratio = "{:10.4f}".format(c.imaginary_light / c.imaginary_heavy) if c.imaginary_light else "{:>10}".format("")
        log.Write(
            "\n    " + c.name[:42].ljust(42) + " " + energy
            + " {:5.0f} {:7.2f} {} {:10.6f}".format(c.degeneracy, 100 * c.population, ratio, c.kie)
            + (" {:10.6f}".format(c.kie_tunnel) if is_kie else "")
        )  # fmt: skip
    log.Write(
        "\n\n  dG: kcal/mol above the lowest conformer of the species; pop: share of the light isotopologue%s; "
        "%s: the conformer against the other side's ensemble"
        % (" (for a transition structure, of its rate, with tunnelling)" if is_kie else "", result.kind)
    )

    log.Write("\n\n" + ("  ").ljust(50))
    log.Write(
        " {:>10} {:>10} {:>10} {:>10} {:>21} {:>6}\n".format(
            result.kind, "1D-tunn", "corr-" + result.kind, "lowest", "range (+/-%.1f)" % result.weight_uncertainty,
            "N_eff",
        )
    )  # fmt: skip
    log.Write(DASH_LINE)
    for r in results:
        log.Write(("\n  " + r.kind + " (ensemble) @ " + str(r.temperature) + " K").ljust(50))
        log.Write(
            " {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}-{:<10.6f} {:6.2f}".format(
                r.kie, r.tunnel_corr, r.kie_tunnel, r.kie_tunnel_lowest, r.kie_tunnel_range[0], r.kie_tunnel_range[1],
                r.n_effective,
            )
        )  # fmt: skip
        if r.reference is not None:
            log.Write(("\n  relative to iso @ " + " / ".join(r.reference.labels) + ":").ljust(50))
            log.Write(" {:10.6f} {:>10} {:10.6f}".format(r.kie_relative, "", r.kie_tunnel_relative))
    log.Write("\n" + DASH_LINE + "\n")
    log.Write(
        "\n  lowest: the lowest conformer of every species alone; range: each free energy moved by +/-%.1f kcal/mol "
        "in turn; N_eff: effective number of %s, 1 / sum(pop^2)\n"
        % (result.weight_uncertainty, "transition-structure conformers" if is_kie else "product conformers")
    )


def _part_name(result):
    """A short description of one step or channel: its transition structure(s)."""
    if isinstance(result, EnsembleIsotopeEffect):
        rows = result.rows("transition structure")
        return rows[0].name + (" (+%d conformers)" % (len(rows) - 1) if len(rows) > 1 else "")
    if isinstance(result, SeriesIsotopeEffect):
        return " -> ".join(_part_name(step) for step in result.steps)
    return result.other.name


def _part_labels(result):
    if isinstance(result, SeriesIsotopeEffect):
        return " | ".join(" / ".join(labels) for labels in result.labels)
    if isinstance(result, EnsembleIsotopeEffect):
        return " / ".join(result.labels)
    return " / ".join(result.reactant.heavy.labels + result.other.heavy.labels)


def _combined_header(log, result, how):
    log.Write(
        "\n\n"
        + (SPACE * 17)
        + "  Temp = "
        + str(result.temperature)
        + "K / Vib. scale factor = "
        + str(result.scale_factor)
        + " / conformer weights: "
        + result.weights
        + " / "
        + how
    )
    if result.tunneling != "bell":
        log.Write(" / tunnelling: " + result.tunneling)
    if result.project:
        log.Write(" / external modes projected")


def _combined_lines(log, results, label, extra_name, extra):
    """The KIE lines of a series or channels, one per temperature, with the reference line."""
    result = results[0]
    log.Write("\n\n" + ("  ").ljust(50))
    log.Write(
        " {:>10} {:>10} {:>10} {:>21} {:>8}\n".format(
            "KIE", "1D-tunn", "corr-KIE", "range (+/-%.1f)" % result.weight_uncertainty, extra_name
        )
    )
    log.Write(DASH_LINE)
    for r in results:
        log.Write(("\n  KIE (" + label + ") @ " + str(r.temperature) + " K").ljust(50))
        spread = (
            "{:10.6f}-{:<10.6f}".format(*r.kie_tunnel_range)
            if r.kie_tunnel_range and r.kie_tunnel_range[0] != r.kie_tunnel_range[1]
            else "{:>21}".format("-")
        )
        value = extra(r)
        log.Write(
            " {:10.6f} {:10.6f} {:10.6f} {} {:>8}".format(
                r.kie, r.tunnel_corr, r.kie_tunnel, spread, "%.4g" % value if value is not None else "-"
            )
        )
        if r.reference is not None:
            log.Write(("\n  relative to its reference:").ljust(50))
            log.Write(" {:10.6f} {:>10} {:10.6f}".format(r.kie_relative, "", r.kie_tunnel_relative))
    log.Write("\n" + DASH_LINE + "\n")


def write_series_results(log, results):
    """Write the steps of a series and its combined KIE (one line per temperature) to the log."""
    results = list(results)
    result = results[0]
    how = {"commitment": "steps weighted by a commitment factor", "user": "steps weighted by given free energies"}
    _combined_header(log, result, how.get(result.weight_source, "steps weighted by computed free energies"))
    log.Write("\n\n  Steps in series at %s K:" % result.temperature)
    log.Write("\n  " + "{:<44} {:>8} {:>8} {:>10} {:>10}".format("", "G", "share %", "KIE", "corr-KIE"))
    for n, step in enumerate(result.steps):
        log.Write("\no step %d: %s" % (n + 1, _part_name(step)))
        log.Write(
            "\n    " + ("iso @ " + _part_labels(step))[:42].ljust(42)
            + " {:8.2f} {:8.2f} {:10.6f} {:10.6f}".format(
                result.free_energies[n], 100 * result.shares[n], step.kie, step.kie_tunnel
            )
        )  # fmt: skip
    log.Write(
        "\n\n  G: effective free energy of the step's transition structure above the lowest (kcal/mol, with "
        "tunnelling); share: of 1/k for the light isotopologue, so the highest transition structure counts most"
    )
    _combined_lines(log, results, "series", "C_f", lambda r: r.commitment)
    if len(result.steps) == 2:
        log.Write("\n  C_f = k2 / k-1, the commitment of the intermediate to going on\n")


def write_channel_results(log, results):
    """Write the channels and their combined KIE (one line per temperature) to the log."""
    results = list(results)
    result = results[0]
    how = {"given": "shares given", "barriers": "shares from given barriers"}
    _combined_header(log, result, how.get(result.share_source, "shares from computed free energies"))
    log.Write("\n\n  Parallel channels at %s K:" % result.temperature)
    log.Write("\n  " + "{:<44} {:>8} {:>10} {:>10}".format("", "share %", "KIE", "corr-KIE"))
    for name, channel, share in zip(result.names, result.channels, result.shares):
        log.Write("\no %s: %s" % (name, _part_name(channel)))
        log.Write(
            "\n    " + ("iso @ " + _part_labels(channel))[:42].ljust(42)
            + " {:8.2f} {:10.6f} {:10.6f}".format(100 * share, channel.kie, channel.kie_tunnel)
        )  # fmt: skip
    log.Write("\n\n  share: of the light isotopologue's rate; 1 / KIE = sum(share / KIE)")
    _combined_lines(log, results, "channels", "s", lambda r: r.selectivity)
    if len(result.channels) == 2:
        log.Write("\n  s = share 2 / share 1, the selectivity between the channels\n")


def write_any(log, results):
    """The results table for whichever kind of result this is."""
    kind = type(results[0])
    writer = {
        EnsembleIsotopeEffect: write_ensemble_results,
        SeriesIsotopeEffect: write_series_results,
        ChannelIsotopeEffect: write_channel_results,
    }.get(kind, write_results)
    writer(log, results)


def read_energies(path):
    """Free energies (and degeneracies) of conformers from a table: file, G or dG, and optionally g.

    Comma, semicolon or whitespace separated; lines starting with # and a header line are skipped. G may be
    '-' to leave a conformer's free energy to Kinisot (a degeneracy alone). Returns {file: (G or None, g)}.
    """
    table = {}
    try:
        with open(path, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
    except OSError as err:
        raise KinisotInputError("cannot read the energies table %s: %s" % (path, err.strerror or err)) from None
    first = True
    for number, line in enumerate(lines, 1):
        fields = line.split("#", 1)[0].replace(",", " ").replace(";", " ").split()
        if not fields:
            continue
        try:
            energy = None if fields[1] == "-" else float(fields[1])
            degeneracy = float(fields[2]) if len(fields) > 2 else 1.0
        except (IndexError, ValueError):
            if first:  # a header
                first = False
                continue
            raise KinisotInputError(
                "%s, line %d: expected a file name, a free energy (or -) and optionally a degeneracy" % (path, number)
            ) from None
        first = False
        table[fields[0]] = (energy, degeneracy)
    if not table:
        raise KinisotInputError("the energies table %s has no rows" % path)
    return table


def _species(groups, table, unit):
    """One entry per species for compute_kie: a file, or Conformers when there are several (or a table row)."""
    species = []
    used = set()
    for files in groups:
        rows = []
        for f in files:
            key = f if f in table else os.path.basename(f) if os.path.basename(f) in table else None
            rows.append(table.get(key))
            used.add(key)
        if len(files) == 1 and rows[0] is None:
            species.append(files[0])
            continue
        energies = [r[0] if r else None for r in rows]
        given = [e is not None for e in energies]
        if any(given) and not all(given):
            missing = [f for f, g in zip(files, given) if not g]
            raise KinisotInputError(
                "the energies table gives free energies for some conformers of a species but not for %s: give all "
                "of them or none" % ", ".join(missing)
            )
        species.append(
            Conformers(
                files,
                free_energies=energies if all(given) else None,
                degeneracy=[r[1] if r else 1.0 for r in rows],
                energy_unit=unit,
            )
        )
    return species, used


def write_json(path, result):
    """Write the full result of one run (or a list, one per temperature) as JSON (overwrites ``path``)."""
    payload = [r.to_dict() for r in result] if isinstance(result, list) else result.to_dict()
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def append_csv(path, result, name=None):
    """Append one row per run to a CSV file (header written when the file is new)."""
    row = result.summary_row()
    if name is not None:
        row = dict(isotopologue=name, **row)
    new_file = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        if new_file:
            writer.writeheader()
        writer.writerow(row)


def build_parser():
    parser = ArgumentParser(
        prog="kinisot",
        description="Kinetic (--ts) and equilibrium (--prd) isotope effects from Gaussian, ORCA or ASE "
        "(machine-learned potential) Hessians, using the Bigeleisen-Mayer equation and a Bell tunnelling correction.",
        epilog="Example: kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 5 -t 393 -s 0.961",
    )
    parser.add_argument(
        "--rct",
        dest="rct",
        action="append",
        nargs="+",
        metavar="FILE",
        help="reactant frequency output (Gaussian .log/.out, ORCA .out/.hess, VibrationsData .json, or a geometry "
        "with --calc); repeat for bimolecular reactions. Several files after one --rct are conformers of that "
        "species, weighted by --weights",
    )
    parser.add_argument(
        "--ts",
        dest="ts",
        action="append",
        nargs="+",
        metavar="FILE",
        help="transition structure frequency output (KIE); several files are conformers",
    )
    parser.add_argument(
        "--prd",
        dest="prd",
        action="append",
        nargs="+",
        metavar="FILE",
        help="product frequency output (EQE); several files are conformers",
    )
    parser.add_argument(
        "--iso",
        dest="label",
        action="append",
        metavar="ATOMS",
        help="atom number(s) to replace with the heavy isotope (2H, 13C, 17O), comma separated, e.g. 7,8. "
        "Give one --iso per file in the order of the --rct then --ts/--prd files, or a single --iso when the "
        "atom numbering is the same in all files. Use 0 for a file without substitution.",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        dest="temperature",
        default=None,
        metavar="K",
        help="temperature in Kelvin (default 298.15); a list (273,298,323) or a range (250:350:10) gives one "
        "result line per temperature",
    )
    parser.add_argument(
        "-s",
        "--scale",
        dest="freq_scale_factor",
        type=float,
        default=None,
        help="vibrational scaling factor (default: ZPE factor from the Truhlar database for the detected level "
        "of theory, else 1.0)",
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
        "--scale-type",
        dest="scale_type",
        choices=["zpe", "harm", "fund"],
        default="zpe",
        help="which Truhlar factor to apply when -s is not given: ZPE (default), harmonic or fundamental",
    )
    parser.add_argument(
        "--tunneling",
        "--tunnelling",
        dest="tunneling",
        choices=["bell", "wigner", "skodje", "none"],
        default="bell",
        help="tunnelling correction for a KIE: Bell infinite parabola (default), Wigner, Skodje-Truhlar (needs the "
        "barrier: --barrier, or electronic energies in the files), or none",
    )
    parser.add_argument(
        "--barrier",
        dest="barrier",
        type=float,
        default=None,
        metavar="KCAL",
        help="barrier height in kcal/mol for --tunneling skodje (default: E(TS) - E(reactants) from the files)",
    )
    parser.add_argument(
        "--project",
        dest="project",
        action="store_true",
        default=None,
        help="project translations and rotations out of the Hessian before diagonalizing (uses the geometry) "
        "instead of discarding the 5/6 lowest modes; the default is on for --calc/ASE inputs and off otherwise",
    )
    parser.add_argument("--no-project", dest="project", action="store_false", help=SUPPRESS)
    parser.add_argument(
        "--calc",
        dest="calc",
        metavar="SPEC",
        help="compute the Hessian of geometry inputs (xyz, extxyz, ... anything ASE reads) with this ASE calculator: "
        "emt, mace_mp[:model], mace_off[:model], mace_omol, orb, sevennet, aimnet2, or module.path:callable. The "
        "Hessian is cached next to the geometry as <name>.hessian.json (needs pip install kinisot[ase])",
    )
    parser.add_argument(
        "--delta",
        dest="delta",
        type=float,
        default=0.01,
        metavar="ANGSTROM",
        help="finite-difference step for --calc Hessians (default 0.01)",
    )
    parser.add_argument(
        "--reference",
        dest="reference",
        action="append",
        metavar="ATOMS",
        help="isotope label(s) of a reference isotopologue (same form as --iso); the KIE is also reported divided "
        "by the reference KIE",
    )
    parser.add_argument(
        "--weights",
        dest="weights",
        choices=WEIGHTS,
        default=None,
        help="how conformers are weighted: qrrho (default: quasi-harmonic free energies, Grimme's entropy "
        "interpolation), rrho (harmonic), user (the free energies in --energies), lowest, or equal",
    )
    parser.add_argument(
        "--energies",
        dest="energies",
        metavar="TABLE",
        help="free energies of conformers: one line per file with its name, G (or dG, or - to compute it) and "
        "optionally a degeneracy; used in place of the computed ones",
    )
    parser.add_argument(
        "--energy-unit",
        dest="energy_unit",
        choices=["kcal/mol", "kJ/mol", "hartree"],
        default="kcal/mol",
        help="unit of the free energies in --energies (default kcal/mol)",
    )
    parser.add_argument(
        "--weight-uncertainty",
        dest="weight_uncertainty",
        type=float,
        default=0.5,
        metavar="KCAL",
        help="the ensemble result gives the range of the KIE when each conformer free energy moves by this much "
        "(default 0.5 kcal/mol)",
    )
    parser.add_argument(
        "--job",
        dest="job",
        metavar="FILE",
        help="a JSON job file with the structures, labels and settings, for transition structures in series, "
        "parallel channels or several isotopologues (docs/job_files.md); the other flags then give only the outputs",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="Kinisot_output.dat",
        metavar="FILE",
        help="results file; new results are appended (default Kinisot_output.dat)",
    )
    parser.add_argument("--overwrite", action="store_true", help="start a fresh results file instead of appending")
    parser.add_argument(
        "--json", dest="json_path", metavar="FILE", help="also write the full result of this run as JSON"
    )
    parser.add_argument(
        "--csv",
        dest="csv_path",
        metavar="FILE",
        help="also append one row with the result of this run to a CSV file (header written when the file is new)",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="do not print results to the terminal")
    parser.add_argument("--version", action="version", version="Kinisot " + __version__)
    return parser


def main(argv=None):
    """Command-line entry point. Returns the process exit code."""
    parser = build_parser()
    options = parser.parse_args(argv)
    if options.job:
        return _main_job(parser, options)
    if not options.rct or not options.label:
        parser.error("--rct and --iso are required (or --job)")
    if options.ts is None and options.prd is None:
        parser.error("either --ts (for a KIE) or --prd (for an EQE) is required")
    if options.ts is not None and options.prd is not None:
        parser.error("--ts and --prd cannot be combined: use --ts for a KIE or --prd for an EQE")
    try:
        temperatures = parse_temperatures(options.temperature or "298.15")
    except ValueError as err:
        parser.error(str(err))
    if min(temperatures) <= 0:
        parser.error("the temperature must be positive")
    if options.freq_scale_factor is not None and options.freq_scale_factor <= 0:
        parser.error("the scaling factor must be positive")

    if options.weight_uncertainty < 0:
        parser.error("--weight-uncertainty cannot be negative")

    is_kie = options.ts is not None
    groups = options.rct + (options.ts if is_kie else options.prd)  # one list of conformer files per species
    # if only one set of labels is provided, assume that the atom numbering is the same for rct and ts/prd
    labels = options.label * 2 if len(options.label) == 1 else list(options.label)
    if len(labels) != len(groups):
        parser.error(
            "%d species were given but %d --iso labels: give one --iso per --rct/--ts/--prd (0 for no substitution) "
            "or a single --iso when the atom numbering is the same" % (len(groups), len(options.label))
        )
    reference = None
    if options.reference:
        reference = options.reference * 2 if len(options.reference) == 1 else list(options.reference)
        if len(reference) != len(groups):
            parser.error("%d species were given but %d --reference labels" % (len(groups), len(options.reference)))
    if options.weights == "user" and not options.energies:
        parser.error("--weights user needs --energies")
    try:
        table = read_energies(options.energies) if options.energies else {}
        rct, used = _species(options.rct, table, options.energy_unit)
        other, used_other = _species(options.ts if is_kie else options.prd, table, options.energy_unit)
    except KinisotInputError as err:
        parser.error(str(err))
    unmatched = sorted(set(table) - used - used_other)
    if unmatched:
        parser.error("the energies table lists files that are not among the inputs: %s" % ", ".join(unmatched))

    try:
        log = Logger(options.output, quiet=options.quiet, overwrite=options.overwrite)
    except OSError as err:
        print(
            "\no  ERROR: cannot open the results file %s: %s\n" % (options.output, err.strerror or err), file=sys.stderr
        )
        return 1
    with log:
        log.Write("\n  " + "KINISOT.py v " + __version__ + ": " + time.strftime("%Y-%m-%d %H:%M") + "\n")
        for files, label in zip(groups, labels):
            names = files[0] if len(files) == 1 else "%s (%d conformers)" % (" ".join(files), len(files))
            log.Write("  Species: {} isotopologue: {}\n".format(names, label))
        try:
            results = []
            seen = set()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", KinisotWarning)
                for temperature in temperatures:
                    results.append(
                        compute_kie(
                            rct=rct,
                            ts=other if is_kie else None,
                            prd=None if is_kie else other,
                            iso=labels,
                            temperature=temperature,
                            scale=options.freq_scale_factor,
                            imag_cutoff=options.freq_cutoff,
                            tunneling=options.tunneling,
                            scale_type=options.scale_type,
                            project=options.project,
                            barrier=options.barrier,
                            reference=reference,
                            calculator=options.calc,
                            delta=options.delta,
                            weights=options.weights,
                            weight_uncertainty=options.weight_uncertainty,
                        )
                    )
            for message in results[0].scaling.messages:
                log.Write("\n  " + message)
            for warning in caught:
                if str(warning.message) not in seen:
                    seen.add(str(warning.message))
                    log.Write("\n  WARNING: " + str(warning.message))
            write_any(log, results)
            if options.json_path:
                write_json(options.json_path, results[0] if len(results) == 1 else results)
            if options.csv_path:
                for result in results:
                    append_csv(options.csv_path, result)
        except KinisotError as err:
            log.Writeonlyfile("o  ERROR: " + str(err))
            print("\no  ERROR: " + str(err) + "\n", file=sys.stderr)
            return 1
        if not options.quiet:
            print("\n  Results appended to " + options.output)
    return 0


# the flags a job file replaces (dest names)
JOB_REPLACES = (
    "rct", "ts", "prd", "label", "temperature", "freq_scale_factor", "freq_cutoff", "scale_type", "tunneling",
    "barrier", "project", "calc", "delta", "reference", "weights", "energies", "energy_unit", "weight_uncertainty",
)  # fmt: skip


def _main_job(parser, options):
    """kinisot --job FILE: every isotopologue of the job at every temperature."""
    from .jobs import load_job, run_job

    defaults = vars(parser.parse_args(["--job", options.job]))
    given = [name for name in JOB_REPLACES if getattr(options, name) != defaults[name]]
    if given:
        parser.error(
            "with --job, the job file gives the structures, labels and settings; the command line gives only the "
            "outputs (-o, --overwrite, -q, --json, --csv). Move these to the job file: %s" % ", ".join(given)
        )
    try:
        job = load_job(options.job)
    except KinisotInputError as err:
        parser.error(str(err))
    try:
        log = Logger(options.output, quiet=options.quiet, overwrite=options.overwrite)
    except OSError as err:
        print(
            "\no  ERROR: cannot open the results file %s: %s\n" % (options.output, err.strerror or err), file=sys.stderr
        )
        return 1
    with log:
        log.Write("\n  " + "KINISOT.py v " + __version__ + ": " + time.strftime("%Y-%m-%d %H:%M") + "\n")
        log.Write("  Job: {} ({})\n".format(options.job, job["kind"].replace("_", " ")))
        try:
            by_name = {}
            seen = set()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", KinisotWarning)
                for temperature in job["temperatures"]:
                    for name, result in run_job(job, temperature):
                        by_name.setdefault(name, []).append(result)
            for warning in caught:
                if str(warning.message) not in seen:
                    seen.add(str(warning.message))
                    log.Write("\n  WARNING: " + str(warning.message))
            for name, results in by_name.items():
                log.Write("\n\n  Isotopologue: " + name)
                write_any(log, results)
            if options.json_path:
                payload = [dict(r.to_dict(), isotopologue=name) for name, rs in by_name.items() for r in rs]
                with open(options.json_path, "w") as handle:
                    json.dump(payload[0] if len(payload) == 1 else payload, handle, indent=2)
                    handle.write("\n")
            if options.csv_path:
                for name, results in by_name.items():
                    for result in results:
                        append_csv(options.csv_path, result, name)
        except KinisotError as err:
            log.Writeonlyfile("o  ERROR: " + str(err))
            print("\no  ERROR: " + str(err) + "\n", file=sys.stderr)
            return 1
        if not options.quiet:
            print("\n  Results appended to " + options.output)
    return 0
