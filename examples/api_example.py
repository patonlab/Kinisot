"""Kinisot from Python: a position scan and a temperature scan of the Claisen KIEs.

Run from the repository root:  python examples/api_example.py
"""

import os

from kinisot import compute_kie, parse_gaussian

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "data", "gaussian")
gs = parse_gaussian(os.path.join(DATA, "claisen_gs.out"))  # parse once, reuse for every substitution
ts = parse_gaussian(os.path.join(DATA, "claisen_ts.out"))

print("13C/2H KIEs for the Claisen rearrangement at 393 K (B3LYP/6-31G(d), scaled by 0.961)")
print("%-8s %10s %10s %10s %10s %10s" % ("atoms", "V-ratio", "ZPE", "TRPF", "KIE", "corr-KIE"))
for atoms in ["1", "2", "3", "4", "5", "6", "7,8"]:
    r = compute_kie(rct=gs, ts=ts, iso=atoms, temperature=393.0, scale=0.961)
    print("%-8s %10.4f %10.4f %10.4f %10.4f %10.4f" % (atoms, r.imag_ratio, r.zpe, r.trpf, r.kie, r.kie_tunnel))

print("\nTemperature dependence of the C4 KIE (Bell and Wigner tunnelling corrections)")
print("%8s %10s %10s %10s" % ("T / K", "KIE", "Bell", "Wigner"))
for temperature in range(300, 501, 50):
    bell = compute_kie(rct=gs, ts=ts, iso="4", temperature=temperature, scale=0.961)
    wigner = compute_kie(rct=gs, ts=ts, iso="4", temperature=temperature, scale=0.961, tunneling="wigner")
    print("%8.0f %10.4f %10.4f %10.4f" % (temperature, bell.kie, bell.kie_tunnel, wigner.kie_tunnel))

# Everything the CLI prints is available on the result object
r = compute_kie(rct=gs, ts=ts, iso="4", temperature=393.0, scale=0.961)
print("\nreaction coordinate: %.1fi (light) / %.1fi (heavy) cm-1" % (r.other.light.imaginary, r.other.heavy.imaginary))
print("external modes discarded from the TS:", ["%.1f" % f for f in r.other.light.species[0].discarded])
print("first JSON keys:", list(r.to_dict())[:6])
