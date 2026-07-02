__version__ = "2.1.0"

from .Hess_to_Freq import read_hess, is_linear, substitute_isotopes  # noqa: E402,F401
from .Kinisot import compute_isotope_effect, calc_rpfr, find_scaling_factor  # noqa: E402,F401
