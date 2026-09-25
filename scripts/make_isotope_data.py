"""Regenerate kinisot/isotope_data.py from the periodictable package (AME 2020 masses).

pip install periodictable
python scripts/make_isotope_data.py
"""

import os

import periodictable as pt

RADIO = {
    "H": [3],
    "C": [11, 14],
    "N": [13],
    "O": [15],
    "F": [18],
    "P": [32, 33],
    "S": [35],
    "Cl": [36],
    "I": [125, 131],
}
HEADER = '''"""Isotope masses (amu) generated from the periodictable package, version %s.

Atomic masses are from the AME 2020 atomic mass evaluation (M. Wang, W. J. Huang,
F. G. Kondev, G. Audi, S. Naimi, Chinese Phys. C 45, 030003 (2021)) and the
natural abundances that pick the most abundant isotope from the IUPAC CIAAW
tables, both as distributed (public domain) by periodictable
(https://github.com/pkienzle/periodictable). Regenerate with
``python scripts/make_isotope_data.py``; do not edit by hand.

ISOTOPE_MASSES maps an element symbol to {mass number: mass} for every
naturally occurring isotope plus a few radioisotopes used as labels (3H, 11C,
14C, 13N, 15O, 18F, 32P, 33P, 35S, 36Cl, 125I, 131I). MOST_ABUNDANT gives the
mass number of the most abundant isotope, which is the light isotopologue
Kinisot builds for every atom (the convention Gaussian uses).
"""
'''


def main():
    lines = [HEADER % pt.__version__, "ISOTOPE_MASSES = {"]
    most = []
    for el in pt.elements:
        if el.number < 1 or el.number > 103:
            continue
        natural = [(i.isotope, i.mass, i.abundance) for i in el if i.abundance > 0]
        if not natural:
            continue  # no natural isotope (Tc, Pm, transuranics): not needed for isotope effects
        extra = [(a, el[a].mass, 0.0) for a in RADIO.get(el.symbol, []) if a in el.isotopes]
        isotopes = sorted({a: m for a, m, _ in natural + extra}.items())
        lines.append('    "%s": {%s},' % (el.symbol, ", ".join("%d: %.9f" % (a, m) for a, m in isotopes)))
        most.append((el.symbol, max(natural, key=lambda t: t[2])[0]))
    lines += ["}", "", "MOST_ABUNDANT = {" + ", ".join('"%s": %d' % (s, a) for s, a in most) + "}", ""]
    target = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "kinisot", "isotope_data.py")
    with open(target, "w") as handle:
        handle.write("\n".join(lines))
    print("wrote", os.path.normpath(target), "with", len(most), "elements")


if __name__ == "__main__":
    main()
