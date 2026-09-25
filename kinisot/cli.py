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
from .exceptions import KinisotError, KinisotWarning

__all__ = ["Logger", "build_parser", "write_results", "write_json", "append_csv", "main"]

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


def write_results(log, result):
    """Write the results table for one calculation to the log."""
    is_kie = result.kind == "KIE"
    rct, oth = result.reactant, result.other
    rct_iso, oth_iso = " / ".join(rct.heavy.labels), " / ".join(oth.heavy.labels)

    log.Write(
        "\n\n"
        + (SPACE * 17)
        + "  Temp = "
        + str(result.temperature)
        + "K / Vib. scale factor = "
        + str(result.scale_factor)
    )
    if is_kie and result.tunneling != "bell":
        log.Write(" / tunnelling: " + result.tunneling)
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
    log.Write(("\n  " + result.kind + " @ " + str(result.temperature) + " K").ljust(50))
    if is_kie:
        log.Write(
            "{:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                result.imag_ratio,
                result.zpe,
                result.exc,
                result.trpf,
                result.kie,
                result.tunnel_corr,
                result.kie_tunnel,
            )
        )
    else:
        log.Write(
            "{:21.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                result.zpe, result.exc, result.trpf, result.kie, result.tunnel_corr, result.kie_tunnel
            )
        )
    log.Write("\n" + DASH_LINE + "\n")

    # Which modes went into the partition functions, so that a misassigned external mode is visible
    log.Write("\n  Vibrational modes (scaled, cm-1): kept in the partition function / discarded as external modes\n")
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


def write_json(path, result):
    """Write the full result of one run as a JSON document (overwrites ``path``)."""
    with open(path, "w") as handle:
        json.dump(result.to_dict(), handle, indent=2)
        handle.write("\n")


def append_csv(path, result):
    """Append one row per run to a CSV file (header written when the file is new)."""
    row = result.summary_row()
    new_file = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        if new_file:
            writer.writeheader()
        writer.writerow(row)


def build_parser():
    parser = ArgumentParser(
        prog="kinisot",
        description="Kinetic (--ts) and equilibrium (--prd) isotope effects from Gaussian or ORCA frequency "
        "calculations, using the Bigeleisen-Mayer equation and a Bell tunnelling correction.",
        epilog="Example: kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 5 -t 393 -s 0.961",
    )
    parser.add_argument(
        "--rct",
        dest="rct",
        action="append",
        required=True,
        metavar="FILE",
        help="reactant frequency output (Gaussian .log/.out, or ORCA .out/.hess); repeat for bimolecular reactions",
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
        help="atom number(s) to replace with the heavy isotope (2H, 13C, 17O), comma separated, e.g. 7,8. "
        "Give one --iso per file in the order of the --rct then --ts/--prd files, or a single --iso when the "
        "atom numbering is the same in all files. Use 0 for a file without substitution.",
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
        choices=["bell", "wigner", "none"],
        default="bell",
        help="tunnelling correction for a KIE: Bell infinite parabola (default), Wigner, or none",
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
            "%d files were given but %d --iso labels: give one --iso per file (0 for no substitution) or a "
            "single --iso when the atom numbering is the same" % (len(files), len(options.label))
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
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", KinisotWarning)
                result = compute_kie(
                    rct=options.rct,
                    ts=options.ts,
                    prd=options.prd,
                    iso=labels,
                    temperature=options.temperature,
                    scale=options.freq_scale_factor,
                    imag_cutoff=options.freq_cutoff,
                    tunneling=options.tunneling,
                    scale_type=options.scale_type,
                )
            for message in result.scaling.messages:
                log.Write("\n  " + message)
            for warning in caught:
                log.Write("\n  WARNING: " + str(warning.message))
            write_results(log, result)
            if options.json_path:
                write_json(options.json_path, result)
            if options.csv_path:
                append_csv(options.csv_path, result)
        except KinisotError as err:
            log.Writeonlyfile("o  ERROR: " + str(err))
            print("\no  ERROR: " + str(err) + "\n", file=sys.stderr)
            return 1
        if not options.quiet:
            print("\n  Results appended to " + options.output)
    return 0
