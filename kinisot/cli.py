"""Command-line interface: argument parsing and .dat/terminal output."""

import json
import sys
import time
from argparse import ArgumentParser

from . import __version__
from .scaling import get_frequency_scaling
from .thermo import compute_isotope_effects

# print formatting
space = "   "
dash = "--"
dash_line = space * 17 + " " + dash * 47


# Enables output to terminal and to text file
class Logger:
    def __init__(self, path, quiet=False):
        self.log = open(path, "w")
        self.quiet = quiet

    # Usable as a context manager so the log file is always closed
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.Finalize()
        return False

    # Write a message to the terminal (unless quiet) and to the log
    def Write(self, message):
        if not self.quiet:
            print(message, end="")
        self.log.write(message)

    # Write a message only to the log and not to the terminal
    def Writeonlyfile(self, message):
        self.log.write("\n" + message + "\n")

    # Write a fatal error, finalize and terminate the program
    def Fatal(self, message):
        print(message + "\n")
        self.log.write(message + "\n")
        self.Finalize()
        sys.exit(1)

    # Finalize the log file
    def Finalize(self):
        if not self.log.closed:
            self.log.close()


def parse_isotopologues(label_flags, error):
    """Split per-file --iso values into named isotopologues.

    Each --iso value belongs to one input file and may hold several
    isotopologues separated by ';', each optionally prefixed 'name='.
    Returns a list of (name, [per-file spec]) with consistent ordering.
    """
    per_file = [[item.strip() for item in flag.split(";")] for flag in label_flags]
    counts = {len(items) for items in per_file}
    if len(counts) != 1:
        error(
            "every --iso flag must define the same number of ';'-separated isotopologues"
        )

    isotopologues = []
    for k in range(counts.pop()):
        name = None
        specs = []
        for items in per_file:
            item = items[k]
            if "=" in item:
                given, _, spec = item.partition("=")
                given = given.strip()
                if name is None:
                    name = given
                elif given and given != name:
                    error(
                        "inconsistent isotopologue names across --iso flags: "
                        "{!r} vs {!r}".format(name, given)
                    )
                item = spec.strip()
            specs.append(item)
        if name is None:
            name = " / ".join(specs)
        isotopologues.append((name, specs))
    return isotopologues


def _stem(path):
    return path.split(".")[0]


def write_detail_block(log, rct, other, is_kie, specs, effect):
    """The classic per-isotopologue block of the .dat output."""
    r = effect.rpfrs
    log.Write("\no " + _stem(rct[0]).ljust(47) + "   " + dash * 47)
    log.Write("\no " + _stem(other[0]).ljust(47))
    if is_kie:
        log.Write("{:10.1f}".format(r[2].im_frequency_wn))
    log.Write(
        "\no "
        + (_stem(rct[0]) + ": iso @ " + " / ".join(specs[0 : len(rct)])).ljust(47)
    )
    import math

    log.Write(
        "           {:10.3e} {:10.3e} {:10.3e}".format(
            math.exp(r[0].ZPE - r[1].ZPE),
            math.exp(r[0].EXC - r[1].EXC),
            math.exp(r[2].PF - r[3].PF),
        )
    )
    log.Write(
        "\no "
        + (_stem(other[0]) + ": iso @ " + " / ".join(specs[len(rct) :])).ljust(47)
    )
    if is_kie:
        log.Write(
            "{:10.1f} {:10.3e} {:10.3e} {:10.3e}".format(
                r[3].im_frequency_wn,
                math.exp(r[2].ZPE - r[3].ZPE),
                math.exp(r[2].EXC - r[3].EXC),
                math.exp(r[0].PF - r[1].PF),
            )
        )
    else:
        log.Write(
            "{:21.3e} {:10.3e} {:10.3e}".format(
                math.exp(r[2].ZPE - r[3].ZPE),
                math.exp(r[2].EXC - r[3].EXC),
                math.exp(r[0].PF - r[1].PF),
            )
        )

    log.Write("\n" + dash_line)
    kind = "KIE" if is_kie else "EQE"
    log.Write(("\n  {} @ {} K".format(kind, effect.temperature)).ljust(50))
    if is_kie:
        log.Write(
            "{:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                effect.v_ratio,
                effect.zpe,
                effect.exc,
                effect.trpf,
                effect.kie,
                effect.bell_correction,
                effect.kie_bell,
                effect.wigner_correction,
                effect.kie_wigner,
            )
        )
    else:
        log.Write(
            "{:21.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                effect.zpe,
                effect.exc,
                effect.trpf,
                effect.kie,
                effect.bell_correction,
                effect.kie_bell,
                effect.wigner_correction,
                effect.kie_wigner,
            )
        )
    log.Write("\n" + dash_line + "\n")


def main():
    # Parse Arguments
    parser = ArgumentParser()
    parser.add_argument(
        "-t",
        dest="temperature",
        action="store",
        type=float,
        default=298.15,
        help="temperature in Kelvin (default 298.15K)",
    )
    parser.add_argument(
        "-s",
        dest="freq_scale_factor",
        action="store",
        type=float,
        default=None,
        help="scale factor for vibrations (default: auto-detect from the level of theory, else 1)",
    )
    parser.add_argument(
        "--iso",
        dest="label",
        action="append",
        required=True,
        help="comma-separated atom number(s) to substitute, each with an optional "
        ":isotope suffix (e.g. '5' or '5:18O,7:2D'; available isotopes: 2D, 3T, "
        "13C, 14C, 15N, 17O, 18O; plain numbers get the default heavy isotope). "
        "Use 0 for a file with no substitution. Repeat the flag once per input "
        "file (a single value is applied to all files). Several isotopologues "
        "may be separated with ';' and optionally named, e.g. 'C5=5;C4=4'",
    )
    parser.add_argument(
        "--ref",
        dest="reference",
        action="store",
        default=None,
        help="isotopologue name to which all other isotope effects are referenced "
        "(divided), matching the internal-standard convention of experimental "
        "KIE measurements",
    )
    parser.add_argument(
        "--cutoff",
        dest="freq_cutoff",
        action="store",
        type=float,
        default=50.0,
        help="frequency cutoff (default = 50 cm-1)",
    )
    parser.add_argument(
        "--rct",
        dest="rct",
        action="append",
        required=True,
        help="reactant output file (repeatable)",
    )
    parser.add_argument(
        "--prd",
        dest="prd",
        action="append",
        help="product output file (for EQE calculation)",
    )
    parser.add_argument(
        "--ts", dest="ts", action="append", help="TS output file (for KIE calculation)"
    )
    parser.add_argument(
        "--output",
        dest="output",
        action="store",
        default="Kinisot_output.dat",
        help="path of the .dat output file (default: Kinisot_output.dat)",
    )
    parser.add_argument(
        "--json",
        dest="json_path",
        action="store",
        default=None,
        help="also write results as JSON to this path",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        action="store_true",
        help="suppress terminal output (the .dat file is still written)",
    )
    parser.add_argument(
        "--version", action="version", version="Kinisot v {}".format(__version__)
    )

    options = parser.parse_args()

    # Validate the run before any output file is created
    if options.ts is None and options.prd is None:
        parser.error(
            "Kinisot requires either a TS (--ts, for a KIE) or a product (--prd, for an EQE)"
        )
    if options.ts is not None and options.prd is not None:
        parser.error("give either --ts (KIE) or --prd (EQE), not both")

    is_kie = options.ts is not None
    other = options.ts if is_kie else options.prd
    files = options.rct + other

    # if only one set of labels is provided, assume that the atom numbering is
    # the same for rct and ts or rct and prd
    if len(options.label) == 1:
        options.label = options.label * 2
    if len(files) != len(options.label):
        parser.error(
            "for multiple reactants you need to specify --iso labels for each of the "
            "{} input files".format(len(files))
        )

    isotopologues = parse_isotopologues(options.label, parser.error)
    names = [name for name, _specs in isotopologues]
    if len(set(names)) != len(names):
        parser.error("isotopologue names must be unique: {}".format(names))
    if options.reference is not None and options.reference not in names:
        parser.error(
            "--ref {!r} does not match any isotopologue name {}".format(
                options.reference, names
            )
        )

    with Logger(options.output, quiet=options.quiet) as log:
        # Write an output file
        log.Write(
            "\n  "
            + "KINISOT.py v "
            + __version__
            + ": "
            + time.strftime("%Y-%m-%d %H:%M")
            + "\n"
        )

        for i, species in enumerate(files):
            log.Write(
                "\n  Species: {} isotopologue(s): {}".format(species, options.label[i])
            )
        log.Write("\n")

        # if not specified try to automatically determine the vibrational scaling factor
        if options.freq_scale_factor is None:
            factor, level, ref, warning = get_frequency_scaling(files)
            if warning is not None:
                log.Writeonlyfile("\nWARNING: " + warning)
            if ref is not None:
                log.Write(
                    "\n  Found vibrational scaling factor {} for {} level of theory".format(
                        factor, level
                    )
                )
                log.Write("\n  REF: " + ref)
            else:
                log.Write(
                    "\n  Unable to find vibrational scaling factor for {}; using value of 1.0".format(
                        level
                    )
                )
            options.freq_scale_factor = factor

        log.Write(
            "\n\n"
            + (space * 17)
            + "  Temp = "
            + str(options.temperature)
            + "K / Vib. scale factor = "
            + str(options.freq_scale_factor)
        )
        log.Write(("\n  ").ljust(50))
        log.Write(
            " {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} \n".format(
                "V-ratio",
                "ZPE",
                "EXC",
                "TRPF",
                "KIE",
                "1D-tunn",
                "corr-KIE",
                "Wigner",
                "W-KIE",
            )
        )

        # Here are the ingredients and final predictions of the isotope effects
        try:
            effects = compute_isotope_effects(
                options.rct,
                options.ts,
                options.prd,
                [specs for _name, specs in isotopologues],
                options.temperature,
                options.freq_scale_factor,
                options.freq_cutoff,
            )
        except (ValueError, FileNotFoundError) as e:
            log.Fatal("\no  " + str(e))

        for (name, specs), effect in zip(isotopologues, effects):
            if len(effects) > 1:
                log.Write("\n  Isotopologue: " + name)
            write_detail_block(log, options.rct, other, is_kie, specs, effect)

        # Summary of referenced (relative) isotope effects
        referenced = {}
        if options.reference is not None:
            ref_effect = effects[names.index(options.reference)]
            log.Write(
                "\n  Isotope effects referenced to {} (reference row shows its absolute value)\n".format(
                    options.reference
                )
            )
            log.Write(
                "  {:20} {:>10} {:>10} {:>10}\n".format(
                    "Isotopologue", "KIE", "corr-KIE", "W-KIE"
                )
            )
            for name, effect in zip(names, effects):
                if name == options.reference:
                    row = (effect.kie, effect.kie_bell, effect.kie_wigner)
                else:
                    row = (
                        effect.kie / ref_effect.kie,
                        effect.kie_bell / ref_effect.kie_bell,
                        effect.kie_wigner / ref_effect.kie_wigner,
                    )
                referenced[name] = row
                log.Write("  {:20} {:10.6f} {:10.6f} {:10.6f}\n".format(name, *row))

        # The 1-D tunneling models are least reliable for primary H/D effects
        for name, effect in zip(names, effects):
            if is_kie and (effect.kie > 1.5 or effect.kie < 0.67):
                log.Write(
                    "\n  WARNING: {} is a large isotope effect; 1-D tunneling "
                    "corrections are least reliable for primary H/D KIEs\n".format(name)
                )

        if options.json_path:
            payload = {
                "kinisot_version": __version__,
                "temperature": options.temperature,
                "freq_scale_factor": options.freq_scale_factor,
                "freq_cutoff": options.freq_cutoff,
                "reactant_files": options.rct,
                ("ts_files" if is_kie else "product_files"): other,
                "reference": options.reference,
                "isotopologues": [
                    dict(
                        name=name,
                        specs=specs,
                        **effect.to_dict(),
                        **(
                            {
                                "referenced_kie": referenced[name][0],
                                "referenced_kie_bell": referenced[name][1],
                                "referenced_kie_wigner": referenced[name][2],
                            }
                            if name in referenced and name != options.reference
                            else {}
                        ),
                    )
                    for (name, specs), effect in zip(isotopologues, effects)
                ],
            }
            with open(options.json_path, "w") as f:
                json.dump(payload, f, indent=2)
            log.Write("\n  JSON results written to " + options.json_path + "\n")


if __name__ == "__main__":
    main()
