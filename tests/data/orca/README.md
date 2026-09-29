# ORCA test fixtures

`claisen_gs.hess` / `claisen_ts.hess` are the Gaussian B3LYP/6-31G(d) Hessians
of the Claisen example rewritten in ORCA's `.hess` layout (`$hessian` in
column blocks of five, `$atoms` with standard atomic weights and Bohr
coordinates, `$vibrational_frequencies`), and the matching `.out` files
carry only the ORCA banner and the `!` keyword line. They exercise the
ORCA reader end to end against the Gaussian golden values; they are not
ORCA output.

The other files are real ORCA 6.1.0 output (Bobby Paton, 2026; first kept on
the `feat/goodvibes-integration` branch), each `.out` with the `.hess` that
ORCA wrote alongside it:

- `pentane_TT` / `pentane_GG`: the anti,anti and gauche,gauche conformers
  of n-pentane, r2SCAN-3c optimization and frequencies (`! r2SCAN-3c OPT
  FREQ`). `tests/test_orca.py` computes the conformational EQE for ²H at
  atom 6 (no scaling factor exists for r2SCAN-3c).
- `hat_gs_freq` / `hat_ts_freq`: the reactant and transition structure of
  a hydrogen-atom transfer, frequency-only runs at the converged geometries
  (the `.inp` files): broken-symmetry M06-2X-D3(0)/6-31+G** with
  SMD(dichloromethane), a triplet guess flipped to M_S = 0. H26 moves from
  C13 to C6 in the transition mode (1974.9i cm⁻¹), so the parabolic-barrier
  crossover temperature is 438 K: at 298 K Kinisot refuses the Bell
  correction, and the test checks the primary and secondary ²H and the ³H
  KIEs, the Wigner and Skodje–Truhlar corrections, and the scaling factor
  (0.968) detected from the `!` line.

The `.out` files are kept whole rather than trimmed to the sections Kinisot
reads. A linear molecule is still wanted (implementation plan, Phase 5).
