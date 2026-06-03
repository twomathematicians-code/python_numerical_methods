# Python Numerical Methods Toolkit

A comprehensive collection of **20 Python programs** implementing fundamental numerical methods, designed as a single, self-contained repository for students, engineers, and researchers.

---

## Table of Contents

- [Numerical Methods Included](#numerical-methods-included)
  - [Root Finding](#root-finding)
  - [Linear Systems](#linear-systems)
  - [Interpolation](#interpolation)
  - [Numerical Integration](#numerical-integration)
  - [Numerical Differentiation](#numerical-differentiation)
  - [Ordinary Differential Equations](#ordinary-differential-equations)
  - [Eigenvalue Problems](#eigenvalue-problems)
- [Project Structure](#project-structure)
- [Running the Programs](#running-the-programs)
- [Requirements](#requirements)
- [Customisation](#customisation)
- [License](#license)

---

## Numerical Methods Included

### Root Finding

| # | Program | Method | Order of Convergence |
|---|---------|--------|---------------------|
| 01 | `01_bisection.py` | Bisection Method | Linear |
| 02 | `02_newton_raphson.py` | Newton-Raphson Method | Quadratic |
| 03 | `03_secant.py` | Secant Method | Superlinear |
| 04 | `04_false_position.py` | False Position (Regula Falsi) | Linear (usually faster than bisection) |

### Linear Systems

| # | Program | Method | Complexity |
|---|---------|--------|-----------|
| 05 | `05_gauss_elimination.py` | Gaussian Elimination (partial pivoting) | O(n^3) |
| 06 | `06_lu_decomposition.py` | LU Decomposition (Doolittle) | O(n^3) factorise, O(n^2) solve |
| 07 | `07_jacobi.py` | Jacobi Iterative Method | O(n^2) per iteration |
| 08 | `08_gauss_seidel.py` | Gauss-Seidel Iterative Method | O(n^2) per iteration |

### Interpolation

| # | Program | Method | Notes |
|---|---------|--------|-------|
| 09 | `09_lagrange_interpolation.py` | Lagrange Polynomial Interpolation | Degree n polynomial |
| 10 | `10_newton_interpolation.py` | Newton's Divided Difference Interpolation | Incremental construction |
| 11 | `11_spline_interpolation.py` | Natural Cubic Spline Interpolation | Piecewise smooth fit |

### Numerical Integration

| # | Program | Method | Accuracy |
|---|---------|--------|---------|
| 12 | `12_trapezoidal.py` | Composite Trapezoidal Rule | O(h^2) |
| 13 | `13_simpson.py` | Composite Simpson's 1/3 Rule | O(h^4) |
| 14 | `14_romberg.py` | Romberg Integration (Richardson extrapolation) | Exponential convergence |

### Numerical Differentiation

| # | Program | Method | Accuracy |
|---|---------|--------|---------|
| 15 | `15_finite_difference.py` | Forward & Central Finite Differences | O(h) / O(h^2) |
| 16 | `16_higher_order_difference.py` | 5-Point Stencil (higher-order) | O(h^4) |

### Ordinary Differential Equations

| # | Program | Method | Accuracy |
|---|---------|--------|---------|
| 17 | `17_euler_ode.py` | Euler's Method | O(h) |
| 18 | `18_rk4_ode.py` | 4th-Order Runge-Kutta | O(h^4) |

### Eigenvalue Problems

| # | Program | Method | Notes |
|---|---------|--------|-------|
| 19 | `19_eigen_power.py` | Power Method (dominant eigenvalue) | Iterative |
| 20 | `20_eigen_qr.py` | QR Algorithm (all eigenvalues) | Modified Gram-Schmidt |

---

## Project Structure

```
python-numerical-methods/
├── README.md                   # This file
├── LICENSE                     # MIT License
├── requirements.txt            # Optional dependencies (numpy, matplotlib)
├── .gitignore
├── Makefile                    # Convenience targets
├── utils.py                    # Shared utilities (print helpers, matrix ops)
├── src/
│   ├── 01_bisection.py
│   ├── 02_newton_raphson.py
│   ├── 03_secant.py
│   ├── 04_false_position.py
│   ├── 05_gauss_elimination.py
│   ├── 06_lu_decomposition.py
│   ├── 07_jacobi.py
│   ├── 08_gauss_seidel.py
│   ├── 09_lagrange_interpolation.py
│   ├── 10_newton_interpolation.py
│   ├── 11_spline_interpolation.py
│   ├── 12_trapezoidal.py
│   ├── 13_simpson.py
│   ├── 14_romberg.py
│   ├── 15_finite_difference.py
│   ├── 16_higher_order_difference.py
│   ├── 17_euler_ode.py
│   ├── 18_rk4_ode.py
│   ├── 19_eigen_power.py
│   └── 20_eigen_qr.py
└── tests/
    └── test_all_methods.py      # Automated verification of all methods
```

---

## Running the Programs

### Run a single method

```bash
# From the project root directory
python src/01_bisection.py

# Or using Make
make run METHOD=01_bisection
```

### Run all methods sequentially

```bash
make run-all
```

### Run the test suite

```bash
make test
# or
python -m pytest tests/
# or (pure stdlib)
python tests/test_all_methods.py
```

### Example output

```
$ python src/01_bisection.py

╔══════════════════════════════════════════════════════════════╗
║  Bisection Method — Root Finding                             ║
╚══════════════════════════════════════════════════════════════╝

Enter lower bound (a): 1
Enter upper bound (b): 2
Enter tolerance: 1e-10
Enter maximum iterations: 100

────────────────────────────────────────────────────────────────
 Iter           a               b           x_mid         f(x_mid)
────────────────────────────────────────────────────────────────
    1    1.00000000    2.00000000    1.50000000   -0.12500000
    2    1.50000000    2.00000000    1.75000000    1.10937500
    ...

Converged after 34 iterations.
Root = 1.5213797068
f(root) = -0.0000000003
```

Each program:
- Is **interactive** — prompts you for input values
- Displays **detailed iteration tables** showing convergence
- Shows **error estimates** comparing to exact solutions (where applicable)
- Has consistent **formatted output** via the shared `utils.py`

---

## Requirements

### Core (zero dependencies — Python 3.7+)

All 20 programs use **only the Python standard library** (`math`, `sys`). No installation is needed.

### Optional (for tests and plots)

```bash
pip install -r requirements.txt
```

| Package | Purpose |
|---------|---------|
| `numpy` | Optional — used in test suite for verification |
| `matplotlib` | Optional — used in test suite for comparison plots |
| `pytest` | Optional — alternative test runner |

---

## Customisation

To implement your own function, simply edit the `f()` function at the top of any source file. For example, in `src/01_bisection.py`:

```python
def f(x):
    # Change this to your own function
    return math.cos(x) - x    # root near 0.7391
```

All root-finding, ODE, integration, and differentiation programs follow this pattern — the target function is defined at the top and clearly labelled.

---

## Comparison with the C++ Version

A companion [C++ implementation](https://github.com/your-username/cpp-numerical-methods) is available with identical algorithms, method-by-method, for performance-critical applications or CUDA/OpenMP parallelisation.

---

## License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details. Free to use, modify, and distribute.
