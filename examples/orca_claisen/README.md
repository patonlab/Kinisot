# ORCA input: the Claisen KIE from `.out` + `.hess` files

ORCA keeps the Hessian in `name.hess` next to `name.out`. Give Kinisot
either file; everything else is as for Gaussian. The fixture files in
`tests/data/orca/` are the Claisen Gaussian Hessians rewritten in ORCA's
layout (see the README there), so the numbers are the same as in the
[claisen](../claisen/README.md) example and check the ORCA reader end to
end. Replace them with your own `name.out`/`name.hess` pair.

```
cd tests/data/orca
kinisot --rct claisen_gs.out --ts claisen_ts.hess --iso 5 -t 393 -s 0.961
kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 393          # scaling factor from the ! line
```

**Result lines** (from `expected_output.dat`):

```
                    V-ratio        ZPE        EXC       TRPF        KIE    1D-tunn   corr-KIE
C5  KIE @ 393.0 K  1.000176   0.999077   1.000680   1.001962   1.001895   1.000044   1.001940
C4  KIE @ 393.0 K  1.012716   1.037214   1.001894   0.978989   1.030281   1.003269   1.033649
```

The second run detects B3LYP/6-31G(d) from the `! B3LYP 6-31G(d) Opt Freq`
line and applies the Truhlar ZPE factor 0.977, hence the different C4 value
from the 0.961-scaled one in the Gaussian example.
