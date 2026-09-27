"""Shi epoxidation over all 18 transition structures: per-structure KIEs and the ensemble KIE.

    python benchmarks/shi_epoxidation/ensemble/analyze.py

Reads ../methylstyrene.log, ../ts10.log (TS A) and ts_<label>.log for the
other 17 transition structures of Singleton & Wang (JACS 2005, 127, 6679).
Prints, for each transition structure, the check against the SI (energy,
zero-point energy, one imaginary mode) and its 13C KIEs next to the
paper's own predictions (SI Table 1, absolute KIEs). It then prints the
ensemble KIEs, relative to the meta carbons as measured, under four
weightings.

The ensemble prototypes IMPLEMENTATION_PLAN.md, Phase 10. With one reactant
conformer,

    KIE = rho_R / sum_j y_j rho'_j

- rho_R is the reactant's reduced partition-function ratio.
- rho'_j is transition structure j's, times (nu_H/nu_L)_j (kappa_H/kappa_L)_j.
- y_j is j's share of the light isotopologue's rate,
  proportional to kappa_L,j exp(-G_j/RT).

Positions that are equivalent in the experiment (the ortho and the meta
carbons) average rho over the pair on each side.
"""

import math
import os
import re
import warnings

import numpy as np

from kinisot import compute_kie, load_hessian
from kinisot.thermo import HARTREE_TO_KCAL_PER_MOL, WAVENUMBER_TO_KELVIN

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPERATURE, SCALE = 273.15, 0.9614
LABELS = ["A", "B", "EA", "DA", "EB", "D", "BD", "DD", "AA", "G", "CA", "H", "E", "C", "AD", "CD", "AB", "CB"]
# [reactant atom, transition-structure atom] (TS A numbering, which every file here uses)
SITES = {
    "C-beta": [("2", "8")],
    "C-alpha": [("1", "7")],
    "CH3": [("3", "9")],
    "ipso": [("4", "2")],
    "ortho": [("10", "3"), ("14", "1")],
    "para": [("12", "5")],
}
META = [("11", "4"), ("13", "6")]
# SI page 80, Table 1: predicted absolute 13C KIEs (C-beta, C-alpha, CH3, ipso, ortho, para).
# Its rows 6-13 are the paper's structures 10-17 (A, B, EA, DA, EB, D, BD, DD).
SI_PREDICTED = {
    "A": (1.022, 1.006, 0.998, 1.001, 1.000, 1.000),
    "B": (1.020, 1.009, 0.999, 1.000, 1.000, 1.000),
    "EA": (1.022, 1.005, 0.998, 1.001, 1.000, 1.000),
    "DA": (1.020, 1.005, 0.998, 1.001, 1.000, 1.000),
    "EB": (1.018, 1.008, 0.999, 1.000, 1.000, 1.000),
    "D": (1.029, 1.006, 0.997, 1.002, 1.000, 1.000),
    "BD": (1.029, 1.006, 0.997, 1.002, 1.000, 1.000),
    "DD": (1.027, 1.005, 0.997, 1.001, 1.000, 1.000),
    "AA": (1.025, 1.007, 0.997, 1.001, 1.000, 1.000),
    "G": (1.036, 1.005, 0.997, 1.002, 1.000, 1.000),
    "CA": (1.026, 1.007, 0.997, 1.001, 1.000, 1.000),
    "H": (1.028, 1.003, 0.998, 1.001, 1.000, 1.000),
    "E": (1.031, 1.004, 0.999, 1.001, 1.000, 1.000),
    "C": (1.014, 1.017, 0.999, 0.999, 1.000, 1.000),
    "AD": (1.026, 1.006, 0.997, 1.001, 1.000, 1.000),
    "CD": (1.027, 1.007, 0.997, 1.001, 1.000, 1.000),
    "AB": (1.021, 1.009, 0.998, 1.001, 1.000, 1.000),
    "CB": (1.023, 1.009, 0.998, 1.001, 1.000, 1.000),
}
# SI page 79: B3LYP/6-311+G** single point + 6-31G* zero-point energy, hartree (the paper's eight structures only)
SI_SINGLE_POINT_ZPE = {
    "A": -1343.404004, "B": -1343.402074, "EA": -1343.399708, "DA": -1343.399115,
    "EB": -1343.398448, "D": -1343.399557, "BD": -1343.399198, "DD": -1343.397175,
}  # fmt: skip
EXPERIMENT = {  # Figure 1 of the paper, two experiments, relative to the meta carbons
    "C-beta": (1.022, 1.020), "C-alpha": (1.005, 1.006), "CH3": (1.002, 1.001),
    "ipso": (0.999, 1.001), "ortho": (1.001, 1.001), "para": (0.998, 1.001),
}  # fmt: skip


def path(label):
    """The frequency job of a transition structure (TS A is the benchmark's ts10.log)."""
    return os.path.join(HERE, "..", "ts10.log") if label == "A" else os.path.join(HERE, "ts_%s.log" % label)


def printed(file, pattern):
    """The last number Gaussian printed after ``pattern``."""
    with open(file, encoding="utf-8") as handle:
        values = re.findall(pattern + r"\s*(-?\d+\.\d+)", handle.read())
    return float(values[-1]) if values else None


def ratios(reactant, ts, pairs, tunneling=True):
    """(rho_R, rho'_TS) for one position, averaged over equivalent atom pairs."""
    r_side, t_side = [], []
    for r_atom, t_atom in pairs:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = compute_kie(
                rct=[reactant], ts=[ts], iso=[r_atom, t_atom], temperature=TEMPERATURE, scale=SCALE,
                tunneling="bell" if tunneling else "none", project=True,
            )  # fmt: skip
        r_side.append(result.reactant.rpfr)
        t_side.append(result.other.rpfr / (result.imag_ratio * result.tunnel_corr))
    return float(np.mean(r_side)), float(np.mean(t_side))


def bell_kappa(imaginary_wn):
    """Bell's infinite-parabola tunnelling factor of one isotopologue."""
    u = WAVENUMBER_TO_KELVIN * imaginary_wn / TEMPERATURE
    return (u / 2) / math.sin(u / 2)


def analyze(tunneling=True):
    """Everything the report prints, as a dictionary (used by tests/test_shi_ensemble.py)."""
    from goodvibes.api import compute_thermo

    reactant = load_hessian(os.path.join(HERE, "..", "methylstyrene.log"))
    rows = {}
    for label in LABELS:
        file = path(label)
        ts = load_hessian(file)
        site = {name: ratios(reactant, ts, pairs, tunneling) for name, pairs in SITES.items()}
        meta = ratios(reactant, ts, META, tunneling)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            probe = compute_kie(rct=[reactant], ts=[ts], iso=["2", "8"], temperature=TEMPERATURE, scale=SCALE,
                                project=True)  # fmt: skip
        thermo = compute_thermo(file, temperature=TEMPERATURE, freq_scale_factor=SCALE, zpe_scale_factor=SCALE)
        with open(file, encoding="utf-8") as handle:
            n_imaginary = [int(n) for n in re.findall(r"\*+\s+(\d+) imaginary frequenc", handle.read())]
        rows[label] = {
            "energy": printed(file, r"SCF Done:  E\(RB3LYP\) ="),
            "zpe": printed(file, r"Zero-point correction="),
            "n_imaginary": n_imaginary[-1] if n_imaginary else 0,
            "imaginary": probe.other.light.imaginary,
            "absolute": {name: r / t for name, (r, t) in site.items()},
            "relative": {name: (r / t) / (meta[0] / meta[1]) for name, (r, t) in site.items()},
            "site": site,
            "meta": meta,
            "e_zpe": printed(file, r"SCF Done:  E\(RB3LYP\) =") + thermo.zpe,
            "g_qrrho": thermo.qh_gibbs_free_energy,
        }
    schemes = {
        "E (B3LYP/6-31G*)": {k: r["energy"] for k, r in rows.items()},
        "E + ZPE": {k: r["e_zpe"] for k, r in rows.items()},
        "qRRHO G, 273 K": {k: r["g_qrrho"] for k, r in rows.items()},
        "6-311+G** + ZPE (SI, 8 structures)": dict(SI_SINGLE_POINT_ZPE),
    }
    ensembles = {}
    for name, free_energy in schemes.items():
        labels = list(free_energy)
        g0 = min(free_energy.values())
        weight = [
            (bell_kappa(rows[k]["imaginary"]) if tunneling else 1.0)
            * math.exp(-(free_energy[k] - g0) * HARTREE_TO_KCAL_PER_MOL / (0.0019872043 * TEMPERATURE))
            for k in labels
        ]
        share = dict(zip(labels, np.array(weight) / sum(weight)))
        meta = ensemble_kie(share, {k: rows[k]["meta"] for k in labels})
        ensembles[name] = {
            "share": share,
            "relative": {s: ensemble_kie(share, {k: rows[k]["site"][s] for k in labels}) / meta for s in SITES},
        }
    return rows, ensembles


def ensemble_kie(share, ratios_by_ts):
    """rho_R / sum_j y_j rho'_j; the reactant's ratio is the same in every transition structure's pair."""
    rho_r = next(iter(ratios_by_ts.values()))[0]
    return rho_r / sum(share[k] * ratios_by_ts[k][1] for k in share)


def main():
    rows, ensembles = analyze()
    experiment = {s: float(np.mean(v)) for s, v in EXPERIMENT.items()}
    print("TS    E (hartree)     nu' (cm-1)  |" + "".join("%9s" % s for s in SITES) + "   SI Table 1 (absolute)")
    for label, row in rows.items():
        calc = [row["absolute"][s] for s in SITES]
        print("%-4s %15.8f %8.1f %2s |%s   %s" % (
            label, row["energy"], row["imaginary"], "i" * row["n_imaginary"],
            "".join("%9.4f" % v for v in calc), " ".join("%.3f" % v for v in SI_PREDICTED[label]),
        ))  # fmt: skip
    print("\nEnsemble KIEs relative to the meta carbons (experiment: mean of two runs)")
    print("%-36s" % "weights" + "".join("%9s" % s for s in SITES) + "   MAD   shares above 5%")
    ensembles["TS A (10) only"] = {"share": {"A": 1.0}, "relative": rows["A"]["relative"]}
    for name, ens in ensembles.items():
        values = [ens["relative"][s] for s in SITES]
        mad = np.mean([abs(v - experiment[s]) for v, s in zip(values, SITES)])
        largest = sorted(ens["share"].items(), key=lambda item: -item[1])
        shares = ", ".join("%s %.0f%%" % (k, 100 * y) for k, y in largest if y > 0.05)
        print("%-36s%s  %.4f  %s" % (name, "".join("%9.4f" % v for v in values), mad, shares))
    print("%-36s%s" % ("experiment", "".join("%9.4f" % experiment[s] for s in SITES)))


if __name__ == "__main__":
    main()
