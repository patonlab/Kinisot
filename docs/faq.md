# FAQ

**Why do my frequencies differ slightly from the ones Gaussian prints?**
Gaussian projects translations and rotations out of the Hessian before
diagonalizing; Kinisot diagonalizes the raw mass-weighted Hessian and drops
the six (five) lowest modes. On a converged geometry the vibrational
frequencies agree to better than 0.05 cm⁻¹ (this is tested), and the
discarded modes are the values on Gaussian's "Low frequencies" line.

**Why not just run Gaussian with `freq=(readisotopes)` for each
isotopologue?** You can, and the KIEs agree to the last printed digit once
the same scaling is used. Kinisot avoids re-running the frequency job for
every label, which matters when scanning all positions of a molecule at
several temperatures.

**What does `--iso 0` mean?** No substituted atom in that file. It is
needed when a side of the reaction is given as several files and only one
of them carries the label. `0` cannot be combined with atom numbers.

**The atom numbering differs between my reactant and my TS.** Give one
`--iso` per file in the order of the `--rct` flags followed by `--ts` or
`--prd`. Kinisot checks that both sides substitute the same elements.

**My TS has two imaginary frequencies.** Kinisot warns and treats the
larger one as the reaction coordinate; the other is discarded with the
external modes, which leaves one rotational residual in the vibrational
product. Re-optimize the transition structure. If the second mode is a
genuine low-frequency torsion below the cutoff, nothing is wrong.

**My reactant has a small imaginary frequency.** Kinisot stops: reactants
and products must be minima. If the mode is a numerical artefact (a few
cm⁻¹ on a floppy molecule), raise `--imag-cutoff` above its magnitude;
otherwise re-optimize.

**Which scaling factor is applied?** The ZPE factor from the Truhlar
database for the detected level of theory, or 1.0 if it is not listed.
Give `-s` to override. Only ZPE, EXC and the tunnelling correction depend
on it; the V-ratio and TRPF are ratios of frequencies and cancel it.

**Why ¹⁷O rather than ¹⁸O?** Historical: the oxygen substitution has
always been ¹⁷O in Kinisot (the isotope measured by ¹⁷O NMR in the Claisen
example). Phase 7 of the implementation plan adds an explicit isotope
syntax (`--iso 5:18O`) and changes the bare-index default to ¹⁸O with a
warning.

**Can I use several conformers?** Not in one run. Compute the KIE for
each reactant/TS conformer pair and Boltzmann-average the rate constants,
or use the lowest-energy pair.

**Does the temperature have to match the Gaussian job?** No. Kinisot
evaluates the partition functions at `-t`; the Hessian does not depend on
temperature.

**How do I get the numbers into a script?** `python -m kinisot ... -q
-o results.dat` and parse the `KIE @` lines, or call
`kinisot.compute_isotope_effect()` from Python. A structured result object
and `--json` output are planned (Phase 4).
