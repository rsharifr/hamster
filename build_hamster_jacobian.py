"""
Build script: read symbolic Jacobian from jacobians/symbolic_hamster.py,
then generate an efficient numerical function in jacobians/hamster_jacobian.py.

Run this after editing symbolic_hamster.py (e.g. when robot physics change).
"""

import os
import sys

# Ensure package is importable when running from project root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from jacobians.symbolic_hamster import J_sym, _output_shape
from jacobians.codegen import matrix_to_flat_code

# Generate efficient numerical code from symbolic J
J_flat_code = matrix_to_flat_code(J_sym, order="F")

lines = [
    '"""',
    "Auto-generated from jacobians/symbolic_hamster.py. Do not edit by hand.",
    "Run build_hamster_jacobian.py to regenerate.",
    '"""',
    "",
    "import numpy as np",
    "",
    "",
    "def hamster_jacobian(D: float, r: float, theta: float) -> np.ndarray:",
    '    """',
    "    Forward Jacobian: body velocity [v_x, v_y, Omega] -> wheel rates [omega1, omega2, omega3].",
    "    Parameters: D = distance center to wheel, r = wheel radius, theta = robot orientation (rad).",
    "    Returns (3, 3) array.",
    '    """',
    "    J_flat = np.array([",
]
for code in J_flat_code:
    lines.append(f"        {code},")
lines.extend([
    "    ], dtype=float)",
    f"    return J_flat.reshape({_output_shape}, order=\"F\")",
    "",
])

out_dir = os.path.join(os.path.dirname(__file__), "jacobians")
out_path = os.path.join(out_dir, "hamster_jacobian.py")
with open(out_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

print(f"Wrote {out_path}")
