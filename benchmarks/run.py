"""Compute every benchmark case and write REPORT.md (and optionally report.json).

    python benchmarks/run.py [--json] [--cases claisen diels_alder]

Each case directory holds a case.json (format in benchmarks/README.md).
Rows whose experimental value is null are computed and listed, but do not
enter the deviation statistics. Cases without structures yet (null
"reactants", or null "transition_structure" and "product"), and entries
without an "iso" label, are listed with their experimental values only.
"""

import argparse
import glob
import json
import math
import os
import sys

from kinisot import __version__, compute_kie

HERE = os.path.dirname(os.path.abspath(__file__))


def load_cases(names=None):
    """Every */case.json (or only the directories in ``names``), with its name and directory attached."""
    cases = []
    for path in sorted(glob.glob(os.path.join(HERE, "*", "case.json"))):
        name = os.path.basename(os.path.dirname(path))
        if names and name not in names:
            continue
        with open(path, encoding="utf-8") as handle:
            case = json.load(handle)
        case["_name"], case["_dir"] = name, os.path.dirname(path)
        cases.append(case)
    return cases


def resolve(case, key):
    """The files listed under ``key`` as paths relative to the case directory, or None."""
    files = case.get(key)
    if not files:
        return None
    return [os.path.normpath(os.path.join(case["_dir"], f)) for f in files]


def measured(entry):
    """Mean of the experimental value(s): a number, or a list of independent measurements."""
    value = entry.get("experimental")
    if isinstance(value, list):
        return sum(value) / len(value)
    return value


def compute(case, iso, reference=None):
    """compute_kie for one isotope label with the case's structures, temperature, scaling and tunnelling."""
    return compute_kie(
        rct=resolve(case, "reactants"),
        ts=resolve(case, "transition_structure"),
        prd=resolve(case, "product"),
        iso=iso,
        temperature=case["temperature"],
        scale=case.get("scale"),
        tunneling=case.get("tunneling", "bell"),
        project=case.get("project"),
        reference=reference,
        imag_cutoff=case.get("imag_cutoff", 50.0),
    )


def geometric_mean(values):
    """Geometric mean, for positions that are equivalent in the experiment."""
    return math.exp(sum(math.log(v) for v in values) / len(values))


def has_structures(case):
    """Reactants and a transition structure (or product) are both given: the case can be computed."""
    return bool(case.get("reactants") and (case.get("transition_structure") or case.get("product")))


def run_case(case):
    """One row per ``kies`` entry: semiclassical and tunnelling-corrected KIE, measurement and deviation."""
    rows = []
    for entry in case["kies"]:
        semiclassical = computed = tunneling = None
        # an entry without an isotopologue label yet is listed with its measurement only
        computable = has_structures(case) and (entry.get("iso") is not None or bool(entry.get("iso_average")))
        if computable and (entry.get("iso_average") or entry.get("reference_average")):
            # positions that are equivalent in the experiment (a rotating methyl group, the two ortho or
            # meta carbons of a phenyl ring) but not in the static structures: geometric mean over them
            results = [compute(case, label) for label in entry.get("iso_average") or [entry["iso"]]]
            computed = geometric_mean([r.kie_tunnel for r in results])
            semiclassical = geometric_mean([r.kie for r in results])
            references = [compute(case, label) for label in entry.get("reference_average", [])]
            if references:
                computed /= geometric_mean([r.kie_tunnel for r in references])
                semiclassical /= geometric_mean([r.kie for r in references])
            tunneling = results[0].tunneling
        elif computable:
            result = compute(case, entry["iso"], entry.get("reference", case.get("reference_isotopologue")))
            computed = result.kie_tunnel_relative if result.reference is not None else result.kie_tunnel
            semiclassical = result.kie_relative if result.reference is not None else result.kie
            tunneling = result.tunneling
        experimental = measured(entry)
        rows.append(
            {
                "position": entry["position"],
                "iso": entry.get("iso"),
                "semiclassical": semiclassical,
                "computed": computed,
                "tunneling": tunneling,
                "experimental": entry.get("experimental"),
                "uncertainty": entry.get("uncertainty"),
                "deviation": (computed - experimental) if None not in (experimental, computed) else None,
                "note": entry.get("note", ""),
            }
        )
    return rows


def format_measurement(value, uncertainty):
    """1.046 ± 0.005, or several independent measurements separated by commas."""
    if not isinstance(value, list):
        value, uncertainty = [value], [uncertainty]
    uncertainty = uncertainty if isinstance(uncertainty, list) else [uncertainty] * len(value)

    def places(x):
        return len(repr(x).partition(".")[2]) if x is not None else 0

    # as many decimals as the value or its uncertainty was given with (JSON drops trailing zeros)
    return ", ".join(
        "%.*f" % (max(places(v), places(u)), v) + (" ± %.*f" % (places(u), u) if u is not None else "")
        for v, u in zip(value, uncertainty)
    )


def format_case(case, rows):
    """The case's section of REPORT.md: sources, conditions, notes and the table of rows."""
    references = case["reference"] if isinstance(case["reference"], list) else [case["reference"]]
    source = " ".join(
        "%s (doi:[%s](https://doi.org/%s)); %s." % (r["citation"], r["doi"], r["doi"], r.get("method", ""))
        for r in references
    )
    if has_structures(case):
        source += " Computed at %s, %s K, scale %s, tunnelling %s%s." % (
            case.get("level_of_theory", "?"),
            case["temperature"],
            case.get("scale") if case.get("scale") is not None else "none",
            case.get("tunneling", "bell"),
            ", relative to isotopologue %s" % case["reference_isotopologue"]
            if case.get("reference_isotopologue")
            else "",
        )
    else:
        source += " Not computed yet."
    lines = ["## %s" % case["name"], "", source, ""]
    if case.get("notes"):
        lines += [case["notes"], ""]
    lines += [
        "| Position | Semiclassical | With tunnelling | Experimental | Deviation | Note |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        exp = (
            "no experimental value"
            if row["experimental"] is None
            else format_measurement(row["experimental"], row["uncertainty"])
        )
        dev = "%+.4f" % row["deviation"] if row["deviation"] is not None else ""
        calc = ["%.4f" % v if v is not None else "–" for v in (row["semiclassical"], row["computed"])]
        lines.append("| %s | %s | %s | %s | %s | %s |" % (row["position"], calc[0], calc[1], exp, dev, row["note"]))
    compared = [r for r in rows if r["deviation"] is not None]
    if compared:
        mad = sum(abs(r["deviation"]) for r in compared) / len(compared)
        lines += ["", "Mean absolute deviation over %d measured positions: %.4f" % (len(compared), mad)]
        if case.get("alternative_to"):
            lines += ["", "An alternative to `%s`, left out of the overall mean." % case["alternative_to"]]
    elif all(r["computed"] is None for r in rows):
        lines += ["", "No structures yet: add the frequency calculations and their paths to `case.json`."]
    else:
        lines += ["", "No experimental values entered yet: enter them in `case.json` from the paper."]
    return "\n".join(lines) + "\n"


def main(argv=None):
    """Run the cases and write REPORT.md (and report.json with --json)."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json", action="store_true", help="also write report.json")
    parser.add_argument("--cases", nargs="*", help="case directory names to run (default: all)")
    options = parser.parse_args(argv)
    cases = load_cases(options.cases)
    if not cases:
        print("no cases found", file=sys.stderr)
        return 1
    report = [
        "# Kinisot benchmark report",
        "",
        "Kinisot %s. Generated by `python benchmarks/run.py`." % __version__,
        "",
    ]
    payload = {}
    all_dev = []
    for case in cases:
        rows = run_case(case)
        payload[case["_name"]] = {"case": {k: v for k, v in case.items() if not k.startswith("_")}, "rows": rows}
        report.append(format_case(case, rows))
        if not case.get("alternative_to"):  # a rejected mechanism, computed for comparison only
            all_dev += [abs(r["deviation"]) for r in rows if r["deviation"] is not None]
    if all_dev:
        report.append(
            "Overall mean absolute deviation over %d measured positions: %.4f\n"
            % (len(all_dev), sum(all_dev) / len(all_dev))
        )
    with open(os.path.join(HERE, "REPORT.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(report))
    if options.json:
        with open(os.path.join(HERE, "report.json"), "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
    print("wrote benchmarks/REPORT.md (%d cases, %d measured positions)" % (len(cases), len(all_dev)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
