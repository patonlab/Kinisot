"""Readers that turn program output into :class:`~kinisot.hessian.HessianInput`.

``load_hessian`` is the single entry point used by the API: it accepts a
ready-made HessianInput or a path and picks the parser from the file
(Gaussian output, ORCA output or ``.hess`` file). ASE follows in Phase 8.
"""

import os

from ..exceptions import KinisotInputError, KinisotParseError
from ..hessian import HessianInput
from .gaussian import parse_gaussian
from .orca import parse_orca

__all__ = ["load_hessian", "detect_program", "parse_gaussian", "parse_orca"]


def detect_program(file):
    """'Gaussian', 'Orca' or 'unknown' from the file name and its first lines."""
    if os.path.splitext(file)[1] == ".hess":
        return "Orca"
    try:
        with open(file, encoding="utf-8", errors="replace") as handle:
            head = [handle.readline() for _ in range(120)]
    except OSError as err:
        raise KinisotParseError("cannot read %s: %s" % (file, err.strerror or err)) from None
    for line in head:
        if "Gaussian" in line:
            return "Gaussian"
        if "* O   R   C   A *" in line:
            return "Orca"
    return "unknown"


def load_hessian(source):
    """Return a HessianInput for ``source`` (a HessianInput or a file path)."""
    if isinstance(source, HessianInput):
        return source
    if isinstance(source, (str, os.PathLike)):
        file = os.fspath(source)
        program = detect_program(file)
        if program == "Orca":
            return parse_orca(file)
        if program == "Gaussian":
            return parse_gaussian(file)
        raise KinisotParseError(
            "%s: not recognized as a Gaussian or ORCA output (ORCA jobs may also be given as the .hess file)" % file
        )
    raise KinisotInputError("cannot read a Hessian from %r: give a file path or a HessianInput" % (source,))
