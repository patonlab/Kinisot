# ORCA test fixtures

`claisen_gs.hess` / `claisen_ts.hess` are the Gaussian B3LYP/6-31G(d) Hessians
of the Claisen example rewritten in ORCA's `.hess` layout (`$hessian` in
column blocks of five, `$atoms` with standard atomic weights and Bohr
coordinates, `$vibrational_frequencies`), and the matching `.out` files
carry only the ORCA banner and the `!` keyword line. They exercise the
ORCA reader end to end against the Gaussian golden values; they are not
ORCA output. Real ORCA fixtures (a small reactant/TS pair and a linear
molecule) are wanted, see the implementation plan, Phase 5.
