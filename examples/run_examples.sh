#!/usr/bin/env bash
# Runs every worked example against the Gaussian outputs in tests/data/gaussian
# and (re)writes examples/<case>/expected_output.dat.
#
#   examples/run_examples.sh            # uses the `kinisot` console script
#   KINISOT="python -m kinisot" examples/run_examples.sh
#
# tests/test_examples.py replays the same commands and compares the result
# lines with the committed expected_output.dat files.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DATA="$HERE/../tests/data/gaussian"
KINISOT="${KINISOT:-kinisot}"

run() {  # run <example directory> <kinisot arguments...>
    local dir="$1"; shift
    $KINISOT "$@" --quiet --output "$HERE/$dir/expected_output.dat"
}

rm -f "$HERE"/*/expected_output.dat
cd "$DATA"

# Claisen rearrangement of allyl vinyl ether at 393 K, B3LYP/6-31G(d) frequencies scaled by 0.961.
# One run per position: 13C at C1-C6, 18O and 17O at O3, and the 2H2 effect at hydrogens 7 and 8.
for atoms in 1 2 3:18O 3:17O 4 5 6 7,8; do
    run claisen --rct claisen_gs.out --ts claisen_ts.out --iso "$atoms" -t 393 -s 0.961
done
# The same C4 KIE with the external modes projected out, relative to C5 as an internal reference,
# with the Skodje-Truhlar correction (barrier from the energies in the files), and as a temperature scan.
run claisen --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 393 -s 0.961 --project
run claisen --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 393 -s 0.961 --reference 5
run claisen --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 393 -s 0.961 --tunneling skodje
run claisen --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 300:400:50 -s 0.961

# Diels-Alder reaction of isoprene with maleic anhydride at 298.15 K, scaled by 0.963.
# First the same 13C KIE from (a) the pre-reaction complex as a single reactant file and
# (b) the two separate reactant files (atom 15 of the TS is atom 6 of the diene).
run diels_alder --rct DATS_rct.out --ts DATS.out --iso 15 -s 0.963
run diels_alder --rct dienophile.out --rct diene.out --ts DATS.out --iso 0 --iso 6 --iso 15 -s 0.963
# Then the remaining diene and dienophile positions (label order: dienophile, diene, TS).
run diels_alder --rct dienophile.out --rct diene.out --ts DATS.out --iso 0 --iso 10 --iso 19 -s 0.963
run diels_alder --rct dienophile.out --rct diene.out --ts DATS.out --iso 0 --iso 1 --iso 10 -s 0.963
run diels_alder --rct dienophile.out --rct diene.out --ts DATS.out --iso 0 --iso 2 --iso 11 -s 0.963
run diels_alder --rct dienophile.out --rct diene.out --ts DATS.out --iso 0 --iso 4 --iso 13 -s 0.963
run diels_alder --rct dienophile.out --rct diene.out --ts DATS.out --iso 1 --iso 0 --iso 1 -s 0.963
run diels_alder --rct dienophile.out --rct diene.out --ts DATS.out --iso 1 --iso 0 --iso 2 -s 0.963
run diels_alder --rct dienophile.out --rct diene.out --ts DATS.out --iso 5 --iso 0 --iso 7 -s 0.963
run diels_alder --rct dienophile.out --rct diene.out --ts DATS.out --iso 5 --iso 0 --iso 5 -s 0.963

# The same Claisen KIE from ORCA-layout files (tests/data/orca, see its README); the second run
# takes the scaling factor from the level of theory on ORCA's ! line.
cd "$HERE/../tests/data/orca"
run orca_claisen --rct claisen_gs.out --ts claisen_ts.hess --iso 5 -t 393 -s 0.961
run orca_claisen --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 393
cd "$DATA"

# The Claisen KIE from Hessians in ASE VibrationsData JSON form (tests/data/ase), without and with projection.
cd "$HERE/../tests/data/ase"
run mlip_claisen --rct claisen_gs.hessian.json --ts claisen_ts.hessian.json --iso 4 -t 393 -s 0.961
run mlip_claisen --rct claisen_gs.hessian.json --ts claisen_ts.hessian.json --iso 4 -t 393 -s 0.961 --project
# The same reaction re-optimized and re-computed with GFN2-xTB (tests/data/xtb, made by scripts/make_claisen_structures.py);
# no scaling factor, external modes projected (the default for Hessians from ASE calculators).
cd "$HERE/../tests/data/xtb"
for atoms in 1 4 6 7,8; do
    run mlip_claisen --rct claisen_gs.hessian.json --ts claisen_ts.hessian.json --iso "$atoms" -t 393
done
cd "$DATA"

# Conformer ensembles (examples/conformers, made by scripts/make_conformer_example.py): GFN2-xTB Hessians of eight
# allyl vinyl ether conformers and the chair and boat Claisen transition structures, at 393 K relative to C5.
# The C4-d2 KIE from the lowest pair, then from every conformer with the degeneracies in degeneracies.txt, then
# every position from the job file.
cd "$HERE/conformers"
run conformers --rct gs_1.hessian.json --ts ts_chair.hessian.json --iso 10,11 -t 393 -s 1 --reference 5
run conformers --rct gs_1.hessian.json gs_2.hessian.json gs_3.hessian.json gs_4.hessian.json gs_5.hessian.json gs_6.hessian.json gs_7.hessian.json gs_8.hessian.json --ts ts_chair.hessian.json ts_boat.hessian.json --iso 10,11 -t 393 -s 1 --reference 5 --energies degeneracies.txt
run conformers --job ensemble.json
cd "$DATA"

# Equilibrium isotope effect: CD3 axial versus equatorial in 1,1,3,3-tetramethylcyclohexane,
# from one frequency calculation with the deuteriums placed on either methyl group. Unscaled.
run eqe_cyclohexane --rct tetramethylcyclohexane.out --prd tetramethylcyclohexane.out --iso 24,25,26 --iso 28,29,30 -s 1 -t 290
run eqe_cyclohexane --rct tetramethylcyclohexane.out --prd tetramethylcyclohexane.out --iso 24,25,26 --iso 28,29,30 -s 1 -t 300

grep -h -E "^  (KIE|EQE)( \((ensemble|series|channels)\))? @|Species:|Isotopologue:" "$HERE"/*/expected_output.dat
