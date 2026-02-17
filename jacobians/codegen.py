"""
Generate efficient numerical Python code from SymPy expressions.

Used by build_hamster_jacobian.py and build_fsquirrle_jacobian.py to produce
the generated jacobian modules from symbolic_hamster.py and symbolic_fsquirrle.py.
"""

import sympy as sp
from sympy.printing.numpy import NumPyPrinter

_printer = NumPyPrinter()


def expr_to_numpy_code(expr: sp.Expr) -> str:
    """Convert a SymPy expression to NumPy code (using np. prefix)."""
    code = _printer.doprint(expr)
    return code.replace("numpy.", "np.")


def matrix_to_flat_code(matrix: sp.Matrix, order: str = "F") -> list[str]:
    """
    Return a list of code strings for each element in column-major (Fortran) order.
    """
    rows, cols = matrix.rows, matrix.cols
    if order == "F":
        indices = [(i, j) for j in range(cols) for i in range(rows)]
    else:
        indices = [(i, j) for i in range(rows) for j in range(cols)]
    return [expr_to_numpy_code(matrix[i, j]) for i, j in indices]


def exprs_to_numpy_code(exprs: list) -> list[str]:
    """Convert a list of SymPy expressions to NumPy code strings."""
    return [expr_to_numpy_code(e) for e in exprs]
