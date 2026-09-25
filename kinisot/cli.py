"""Console entry point (``kinisot`` command).

The argument parsing and formatting still live in :mod:`kinisot.Kinisot`;
Phase 4 of the implementation plan moves them here.
"""

from .Kinisot import build_parser, main

__all__ = ["build_parser", "main"]
