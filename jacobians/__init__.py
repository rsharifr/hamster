"""
Jacobians package: numerical Jacobian and basis functions.

Flow:
  1. Edit symbolic kinematics (physics) in:
     - jacobians/symbolic_hamster.py
     - jacobians/symbolic_fsquirrle.py
  2. Regenerate efficient numerical code:
     - python build_hamster_jacobian.py
     - python build_fsquirrle_jacobian.py
  3. Main scripts (hamster.py, fsquirrle.py) import from here.

This package does not require sympy at runtime; only the build scripts do.
"""

from .hamster_jacobian import hamster_jacobian
from .fsquirrle_jacobian import (
    fsquirrle_jacobian_inverse,
    fsquirrle_jacobian,
    fsquirrle_bases,
)

__all__ = [
    "hamster_jacobian",
    "fsquirrle_jacobian_inverse",
    "fsquirrle_jacobian",
    "fsquirrle_bases",
]
