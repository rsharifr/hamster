## Python port of HAMSTER MATLAB toolbox

This folder contains a Python reimplementation of the original HAMSTER MATLAB
codebase. It mirrors the functionality of the MATLAB scripts:

- `run__HAMSTER_motorSelection.m` → `hamster.run_hamster_motor_selection`
- `run__FSquirrle_motorSelection.m` → `fsquirrle.run_fsquirrle_motor_selection`
- `run__1D_lqr_test.m` → `lqr_test.run_1d_lqr_test`

### Installation

From this folder (or repo root), create a virtual environment and install
dependencies:

```bash
python -m venv .venv
.venv\Scripts\activate  # on Windows
pip install -r requirements.txt
```

### Modules

- `hamster.py`
  - `hamster_jacobian(D, r, theta)` – Jacobian used by the HAMSTER motor selection.
  - `run_hamster_motor_selection(...)` – replicates the torque/speed sweeps and plots.

- `fsquirrle.py`
  - `fibonacci_sphere(n)` – direction sampling on the unit sphere.
  - `fsquirrle_jacobian_inverse(...)` / `fsquirrle_jacobian(...)` – kinematic Jacobians.
  - `fsquirrle_bases(...)` – cable force bases.
  - `run_fsquirrle_motor_selection(...)` – full workspace motor selection analysis.

- `lqr_test.py`
  - `lqr(A, B, Q, R)` – simple continuous-time LQR.
  - `run_1d_lqr_test(...)` – 1D admittance/LQR test with noise and delay.

### Examples

Run the Python equivalents from the repo root:

```bash
python -m python_port.hamster
python -m python_port.fsquirrle
python -m python_port.lqr_test
```

Or from a Python REPL:

```python
from python_port import hamster, fsquirrle, lqr_test

results_hamster = hamster.run_hamster_motor_selection()
results_fs = fsquirrle.run_fsquirrle_motor_selection()
results_lqr = lqr_test.run_1d_lqr_test()
```

