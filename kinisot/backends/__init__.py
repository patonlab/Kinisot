"""Readers that turn program output into :class:`~kinisot.hessian.HessianInput`.

``load_hessian`` is the single entry point used by the API: it accepts a
ready-made HessianInput or a path and picks the parser. Only Gaussian is
implemented here; ORCA (via GoodVibes) and ASE follow in later phases.
"""

import os

from ..exceptions import KinisotInputError
from ..hessian import HessianInput
from .gaussian import parse_gaussian

__all__ = ["load_hessian", "parse_gaussian"]


def load_hessian(source):
    """Return a HessianInput for ``source`` (a HessianInput or a file path)."""
    if isinstance(source, HessianInput):
        return source
    if isinstance(source, (str, os.PathLike)):
        return parse_gaussian(os.fspath(source))
    raise KinisotInputError("cannot read a Hessian from %r: give a file path or a HessianInput" % (source,))
