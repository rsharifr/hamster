"""
Python port of the HAMSTER MATLAB toolbox.

Layout:
- build_hamster_jacobian.py / build_fsquirrle_jacobian.py: derive Jacobians
  (run these first to generate jacobians/*.py).
- jacobians/: generated numerical Jacobian and basis functions.
- hamster.py, fsquirrle.py, lqr_test.py: main scripts that use the Jacobians.

Usage: run build_*_jacobian.py to regenerate Jacobians; then run hamster.py,
fsquirrle.py, or lqr_test.py as needed.
"""

