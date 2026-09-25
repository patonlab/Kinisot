# ASE test fixtures

`claisen_gs.hessian.json` / `claisen_ts.hessian.json` are the Gaussian
B3LYP/6-31G(d) Hessians of the Claisen example written in the format Kinisot
uses for ASE inputs (`ase.vibrations.VibrationsData.todict()` encoded with
`ase.io.jsonio`, geometry in Angstrom, Hessian in eV/Angstrom^2, the
electronic energy in eV under `atoms.info["energy"]`). They exercise the
ASE reader and the unit conversions against the Gaussian golden values.
`kinisot.save_hessian_json` writes this format; `--calc` writes it as the
cache next to a geometry file.
