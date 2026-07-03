"""Command-line interface: argument parsing and terminal/.dat output."""

import json
import math
import sys
import time
import warnings
from argparse import ArgumentParser

from rich import box
from rich.console import Console
from rich.table import Table

from . import __version__
from .scaling import get_frequency_scaling
from .thermo import compute_isotope_effects


# Legacy plain-text logger, kept for the kinisot.Kinisot compatibility shim
class Logger:
    def __init__(self, path, quiet=False):
        self.log = open(path, "w")
        self.quiet = quiet

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.Finalize()
        return False

    def Write(self, message):
        if not self.quiet:
            print(message, end="")
        self.log.write(message)

    def Writeonlyfile(self, message):
        self.log.write("\n" + message + "\n")

    def Fatal(self, message):
        print(message + "\n")
        self.log.write(message + "\n")
        self.Finalize()
        sys.exit(1)

    def Finalize(self):
        if not self.log.closed:
            self.log.close()


class RichLogger:
    """Renders every message/table to the terminal and to the .dat file."""

    def __init__(self, path, quiet=False):
        self._handle = open(path, "w")
        self.terminal = Console(quiet=quiet, highlight=False)
        self.file = Console(file=self._handle, width=150, highlight=False)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    def print(self, *args, **kwargs):
        self.terminal.print(*args, **kwargs)
        self.file.print(*args, **kwargs)

    def warn(self, message):
        self.print("[yellow]WARNING:[/yellow] " + message)
        self.print()

    def fatal(self, message):
        # always reaches the terminal, even in quiet mode
        print(message)
        self.file.print(message)
        self.close()
        sys.exit(1)

    def close(self):
        if not self._handle.closed:
            self._handle.close()


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


def _short(path):
    return path.rsplit("/", 1)[-1]


def _num(value, digits=6):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "[dim]—[/dim]"
    return "{:.{}f}".format(value, digits)


def species_table(rct, other, is_kie, effects):
    table = Table(box=box.SIMPLE, title=None, show_edge=False)
    table.add_column("Role", style="bold")
    table.add_column("File")
    table.add_column("ν‡ / cm⁻¹", justify="right")
    light = effects[0].rpfrs
    for f in rct:
        table.add_row("Reactant", _short(f), "")
    for f in other:
        table.add_row(
            "TS" if is_kie else "Product",
            _short(f),
            "{:.1f}i".format(light[2].im_frequency_wn) if is_kie else "",
        )
    return table


def results_table(names, effects, is_kie, temperature, scale):
    kind = "Kinetic" if is_kie else "Equilibrium"
    table = Table(
        title="{} isotope effects @ {} K (vib. scale factor {})".format(
            kind, temperature, scale
        ),
        box=box.SIMPLE_HEAVY,
    )
    table.add_column("Isotopologue", style="bold")
    if is_kie:
        for col in (
            "V-ratio",
            "ZPE",
            "EXC",
            "TRPF",
            "KIE",
            "κ Wigner",
            "KIE×κW",
            "κ Bell",
            "KIE×κB",
        ):
            table.add_column(col, justify="right")
        for name, e in zip(names, effects):
            table.add_row(
                name,
                _num(e.v_ratio),
                _num(e.zpe),
                _num(e.exc),
                _num(e.trpf),
                "[bold]{}[/bold]".format(_num(e.kie)),
                _num(e.wigner_correction),
                _num(e.kie_wigner),
                _num(e.bell_correction),
                _num(e.kie_bell),
            )
    else:
        for col in ("ZPE", "EXC", "TRPF", "EQE"):
            table.add_column(col, justify="right")
        for name, e in zip(names, effects):
            table.add_row(
                name,
                _num(e.zpe),
                _num(e.exc),
                _num(e.trpf),
                "[bold]{}[/bold]".format(_num(e.kie)),
            )
    return table


def components_table(names, effects, is_kie):
    """Per-species Bigeleisen-Mayer ingredients (light/heavy RPFR ratios)."""
    table = Table(
        title="RPFR components (light/heavy ratios per species)", box=box.SIMPLE
    )
    table.add_column("Isotopologue", style="bold")
    table.add_column("Species")
    for col in ("ZPE ratio", "EXC ratio", "PF ratio"):
        table.add_column(col, justify="right")
    for name, e in zip(names, effects):
        r = e.rpfrs
        table.add_row(
            name,
            "reactant(s)",
            "{:.3e}".format(math.exp(r[0].ZPE - r[1].ZPE)),
            "{:.3e}".format(math.exp(r[0].EXC - r[1].EXC)),
            "{:.3e}".format(math.exp(r[1].PF - r[0].PF)),
        )
        table.add_row(
            "",
            "TS" if is_kie else "product(s)",
            "{:.3e}".format(math.exp(r[2].ZPE - r[3].ZPE)),
            "{:.3e}".format(math.exp(r[2].EXC - r[3].EXC)),
            "{:.3e}".format(math.exp(r[3].PF - r[2].PF)),
        )
    return table


def referenced_table(names, effects, reference):
    ref_effect = effects[names.index(reference)]
    table = Table(
        title="Isotope effects referenced to {} "
        "(reference row shows its absolute value)".format(reference),
        box=box.SIMPLE_HEAVY,
    )
    table.add_column("Isotopologue", style="bold")
    for col in ("KIE", "KIE×κW", "KIE×κB"):
        table.add_column(col, justify="right")
    rows = {}
    for name, e in zip(names, effects):
        if name == reference:
            row = (e.kie, e.kie_wigner, e.kie_bell)
        else:
            row = (
                e.kie / ref_effect.kie,
                e.kie_wigner / ref_effect.kie_wigner,
                e.kie_bell / ref_effect.kie_bell,
            )
        rows[name] = row
        table.add_row(name, *[_num(v) for v in row])
    return table, rows


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

    with RichLogger(options.output, quiet=options.quiet) as log:
        log.print(
            "[bold]KINISOT[/bold] v {}  ·  {}".format(
                __version__, time.strftime("%Y-%m-%d %H:%M")
            )
        )
        log.print()

        # if not specified try to automatically determine the vibrational scaling factor
        if options.freq_scale_factor is None:
            factor, level, ref, warning = get_frequency_scaling(files)
            if warning is not None:
                log.warn(warning)
            if ref is not None:
                log.print(
                    "Vibrational scaling factor [bold]{}[/bold] for {} "
                    "level of theory".format(factor, level)
                )
                log.print("[dim]REF: {}[/dim]".format(ref))
            else:
                log.print(
                    "No vibrational scaling factor found for {}; using 1.0".format(
                        level
                    )
                )
            options.freq_scale_factor = factor

        # Compute all isotopologues, collecting library warnings for tidy display
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
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
                log.fatal("ERROR: " + str(e))

        log.print()
        log.print(species_table(options.rct, other, is_kie, effects))
        log.print()
        for message in dict.fromkeys(str(w.message) for w in caught):
            log.warn(message)

        log.print(
            results_table(
                names, effects, is_kie, options.temperature, options.freq_scale_factor
            )
        )
        log.print(components_table(names, effects, is_kie))

        referenced = {}
        if options.reference is not None:
            table, referenced = referenced_table(names, effects, options.reference)
            log.print(table)

        # The 1-D tunneling models are least reliable for primary H/D effects
        for name, effect in zip(names, effects):
            if is_kie and (effect.kie > 1.5 or effect.kie < 0.67):
                log.warn(
                    "{} is a large isotope effect; 1-D tunneling corrections "
                    "are least reliable for primary H/D KIEs".format(name)
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
                                "referenced_kie_bell": referenced[name][2],
                                "referenced_kie_wigner": referenced[name][1],
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
            log.print("[dim]JSON results written to {}[/dim]".format(options.json_path))


if __name__ == "__main__":
    main()
