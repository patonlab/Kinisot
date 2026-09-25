"""Exceptions and warnings raised by Kinisot's library code.

Library functions never call ``sys.exit``; they raise one of these and the
command-line entry point turns them into messages and exit codes. Both
error classes also derive from ``ValueError`` so that code written against
Kinisot 2.0 (``except ValueError``) keeps working.
"""


class KinisotError(Exception):
    """Base class for all Kinisot errors."""


class KinisotParseError(KinisotError, ValueError):
    """An input file could not be read or does not contain what Kinisot needs."""


class KinisotInputError(KinisotError, ValueError):
    """The requested calculation is inconsistent (labels, files, structures)."""


class KinisotWarning(UserWarning):
    """Non-fatal problems that the user should see (e.g. extra imaginary modes)."""
