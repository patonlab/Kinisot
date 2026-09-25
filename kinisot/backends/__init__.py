"""Readers that turn program output into :class:`~kinisot.hessian.HessianInput`.

``load_hessian`` is the single entry point used by the API: it accepts a
ready-made HessianInput or a path and picks the parser from the file
(Gaussian output, ORCA output or ``.hess`` file, a VibrationsData JSON file,
or any geometry ASE can read when a calculator is given).
"""

import os

from ..exceptions import KinisotInputError, KinisotParseError
from ..hessian import HessianInput
from .ase import hessian_for_geometry, is_ase_json, parse_ase_json
from .gaussian import parse_gaussian
from .orca import parse_orca

__all__ = ["load_hessian", "detect_program", "parse_gaussian", "parse_orca", "parse_ase_json"]


def detect_program(file):
    """'Gaussian', 'Orca', 'ase' (VibrationsData JSON) or 'unknown' from the file name and its first lines."""
    if os.path.splitext(file)[1] == ".hess":
        return "Orca"
    if is_ase_json(file):
        return "ase"
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


def load_hessian(source, calculator=None, delta=0.01):
    """Return a HessianInput for ``source``.

    ``source`` is a HessianInput, a Gaussian or ORCA output, an ORCA ``.hess``
    file, a VibrationsData JSON file, or (with ``calculator``, a ``--calc``
    specification or an ASE calculator object) a geometry file whose Hessian
    is computed and cached next to it.
    """
    if isinstance(source, HessianInput):
        return source
    if not isinstance(source, (str, os.PathLike)):
        raise KinisotInputError("cannot read a Hessian from %r: give a file path or a HessianInput" % (source,))
    file = os.fspath(source)
    program = detect_program(file)
    if program == "Orca":
        return parse_orca(file)
    if program == "Gaussian":
        return parse_gaussian(file)
    if program == "ase":
        return parse_ase_json(file)
    if calculator is not None:
        if isinstance(calculator, str):
            return hessian_for_geometry(file, calculator, delta=delta)
        from ase.io import read

        from .ase import hessian_from_calculator

        atoms = read(file)
        return hessian_from_calculator(
            atoms[-1] if isinstance(atoms, list) else atoms, calculator, delta=delta, source=file
        )
    raise KinisotParseError(
        "%s: not recognized as a Gaussian or ORCA output or a VibrationsData JSON file (ORCA jobs may also be given "
        "as the .hess file; a geometry file needs --calc to compute its Hessian)" % file
    )
