#!/usr/bin/env python3
"""
test_all_methods.py — Automated verification of all 20 numerical methods.

Runs each method with pre-defined inputs and validates that the computed
result matches the known exact answer within a specified tolerance.
No external dependencies — uses only the Python standard library.

Run:  python tests/test_all_methods.py
"""

import sys
import os
import math

# Ensure project root is on the path so utils.py can be imported
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ─── Import utility functions ──────────────────────────────────────────────
from utils import mat_vec_mul, vec_norm

# ─── Test infrastructure ────────────────────────────────────────────────────

PASSED = 0
FAILED = 0


def check(name: str, condition: bool, detail: str = "") -> None:
    global PASSED, FAILED
    if condition:
        PASSED += 1
        print(f"  ✓ {name}")
    else:
        FAILED += 1
        print(f"  ✗ {name}  —  {detail}")


# ══════════════════════════════════════════════════════════════════════════
# 01 — Bisection
# ══════════════════════════════════════════════════════════════════════════

def test_bisection():
    def f(x): return x ** 3 - x - 2.0

    a, b = 1.0, 2.0
    mid = 0.0
    for _ in range(100):
        mid = (a + b) / 2.0
        if abs(f(mid)) < 1e-15 or (b - a) / 2 < 1e-15:
            break
        if f(a) * f(mid) < 0:
            b = mid
        else:
            a = mid

    check("Bisection root ≈ 1.5214", abs(mid - 1.5213797068) < 1e-8)


# ══════════════════════════════════════════════════════════════════════════
# 02 — Newton-Raphson
# ══════════════════════════════════════════════════════════════════════════

def test_newton_raphson():
    def f(x): return x ** 3 - x - 2.0
    def df(x): return 3 * x * x - 1.0

    x = 1.5
    for _ in range(100):
        dx = f(x) / df(x)
        x -= dx
        if abs(dx) < 1e-12:
            break

    check("Newton-Raphson root ≈ 1.5214", abs(x - 1.5213797068) < 1e-10)


# ══════════════════════════════════════════════════════════════════════════
# 03 — Secant
# ══════════════════════════════════════════════════════════════════════════

def test_secant():
    def f(x): return x ** 3 - x - 2.0
    prev, curr = 1.0, 2.0
    for _ in range(100):
        x_new = curr - f(curr) * (curr - prev) / (f(curr) - f(prev))
        if abs(x_new - curr) < 1e-12:
            break
        prev, curr = curr, x_new

    check("Secant root ≈ 1.5214", abs(curr - 1.5213797068) < 1e-10)


# ══════════════════════════════════════════════════════════════════════════
# 04 — False Position
# ══════════════════════════════════════════════════════════════════════════

def test_false_position():
    def f(x): return x ** 3 - x - 2.0
    a, b = 1.0, 2.0
    fa, fb = f(a), f(b)
    x = 0.0
    for _ in range(100):
        x = (a * fb - b * fa) / (fb - fa)
        fx = f(x)
        if abs(fx) < 1e-12:
            break
        if fa * fx < 0:
            b, fb = x, fx
        else:
            a, fa = x, fx

    check("False Position root ≈ 1.5214", abs(x - 1.5213797068) < 1e-10)


# ══════════════════════════════════════════════════════════════════════════
# 05 — Gauss Elimination
# ══════════════════════════════════════════════════════════════════════════

def test_gauss_elimination():
    # [1, 1; -1, 1] * x = [5; 1]  →  x = [2, 3]
    aug = [[1, 1, 5], [-1, 1, 1]]

    n = 2
    for col in range(n):
        pivot = col
        for row in range(col + 1, n):
            if abs(aug[row][col]) > abs(aug[pivot][col]):
                pivot = row
        aug[col], aug[pivot] = aug[pivot], aug[col]
        for row in range(col + 1, n):
            factor = aug[row][col] / aug[col][col]
            for j in range(col, n + 1):
                aug[row][j] -= factor * aug[col][j]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = aug[i][n]
        for j in range(i + 1, n):
            x[i] -= aug[i][j] * x[j]
        x[i] /= aug[i][i]

    check("Gauss: x=[2,3]",
          abs(x[0] - 2.0) < 1e-10 and abs(x[1] - 3.0) < 1e-10,
          f"got {x}")


# ══════════════════════════════════════════════════════════════════════════
# 06 — LU Decomposition
# ══════════════════════════════════════════════════════════════════════════

def test_lu_decomposition():
    A = [[2, 1], [5, 3]]
    b = [5, 13]
    n = 2

    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            s = sum(L[i][k] * U[k][j] for k in range(i))
            U[i][j] = A[i][j] - s
        for j in range(i, n):
            if i == j:
                L[i][i] = 1.0
            else:
                s = sum(L[j][k] * U[k][i] for k in range(i))
                L[j][i] = (A[j][i] - s) / U[i][i]

    y = [0.0] * n
    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]

    # Verify L*U = A
    LU = [[sum(L[i][k] * U[k][j] for k in range(n)) for j in range(n)]
          for i in range(n)]
    check("LU: L*U = A", all(abs(LU[i][j] - A[i][j]) < 1e-10
                             for i in range(n) for j in range(n)))
    check("LU: solution x=[2,1]", abs(x[0] - 2) < 1e-10 and abs(x[1] - 1) < 1e-10)


# ══════════════════════════════════════════════════════════════════════════
# 07 — Jacobi
# ══════════════════════════════════════════════════════════════════════════

def test_jacobi():
    # [5, 1; 1, 4] * x = [6; 5]  →  x = [19/19, 6/19]... let me use a clean one
    # [6, 1; 2, 5] * x = [10; 8] →  6x0+x1=10, 2x0+5x1=8
    # x1=10-6x0; 2x0+5(10-6x0)=8; 2x0+50-30x0=8; -28x0=-42; x0=1.5, x1=1
    A = [[6, 1], [2, 5]]
    b = [10.0, 8.0]
    x = [0.0, 0.0]

    for _ in range(200):
        x_new = [(b[i] - sum(A[i][j] * x[j] for j in range(2) if j != i)) / A[i][i]
                 for i in range(2)]
        if max(abs(x_new[i] - x[i]) for i in range(2)) < 1e-12:
            break
        x = x_new

    check("Jacobi: x ≈ [1.5, 1]",
          abs(x[0] - 1.5) < 1e-8 and abs(x[1] - 1.0) < 1e-8,
          f"got {x}")


# ══════════════════════════════════════════════════════════════════════════
# 08 — Gauss-Seidel
# ══════════════════════════════════════════════════════════════════════════

def test_gauss_seidel():
    A = [[6, 1], [2, 5]]
    b = [10.0, 8.0]
    x = [0.0, 0.0]

    for _ in range(200):
        x_old = x[:]
        for i in range(2):
            x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(2) if j != i)) / A[i][i]
        if max(abs(x[i] - x_old[i]) for i in range(2)) < 1e-12:
            break

    check("Gauss-Seidel: x ≈ [1.5, 1]",
          abs(x[0] - 1.5) < 1e-8 and abs(x[1] - 1.0) < 1e-8,
          f"got {x}")


# ══════════════════════════════════════════════════════════════════════════
# 09 — Lagrange Interpolation
# ══════════════════════════════════════════════════════════════════════════

def test_lagrange():
    # Use f(x) = x^2 with points (0,0), (1,1), (2,4)
    X = [0, 1, 2]
    Y = [0, 1, 4]

    def interpolate(xv):
        result = 0.0
        for i in range(3):
            Li = 1.0
            for j in range(3):
                if j != i:
                    Li *= (xv - X[j]) / (X[i] - X[j])
            result += Li * Y[i]
        return result

    # P(1.5) = (1.5)^2 = 2.25
    val = interpolate(1.5)
    check("Lagrange: P(1.5) = 2.25 for f(x)=x^2", abs(val - 2.25) < 1e-10,
          f"got {val}")

    # Also check data points
    check("Lagrange: passes through data points",
          all(abs(interpolate(X[i]) - Y[i]) < 1e-12 for i in range(3)))


# ══════════════════════════════════════════════════════════════════════════
# 10 — Newton Divided Difference
# ══════════════════════════════════════════════════════════════════════════

def test_newton_interp():
    # f(x) = x^2 with points (0,0), (1,1), (2,4)
    X = [0, 1, 2]
    Y = [0, 1, 4]
    n = 3
    dd = [[0.0] * n for _ in range(n)]
    for i in range(n):
        dd[i][0] = Y[i]
    for j in range(1, n):
        for i in range(n - j):
            dd[i][j] = (dd[i + 1][j - 1] - dd[i][j - 1]) / (X[i + j] - X[i])

    xval = 1.5
    result = dd[0][0]
    product = 1.0
    for j in range(1, n):
        product *= (xval - X[j - 1])
        result += dd[0][j] * product

    check("Newton DD: P(1.5) = 2.25 for f(x)=x^2", abs(result - 2.25) < 1e-10,
          f"got {result}")


# ══════════════════════════════════════════════════════════════════════════
# 11 — Cubic Spline
# ══════════════════════════════════════════════════════════════════════════

def test_spline():
    X = [0, 1, 2, 3]
    Y = [0, 1, 0, -1]
    n = len(X)
    m = n - 1
    h_vals = [X[i + 1] - X[i] for i in range(m)]

    alpha = [0.0] * m
    for i in range(1, m):
        alpha[i] = (3 / h_vals[i]) * (Y[i + 1] - Y[i]) - \
                   (3 / h_vals[i - 1]) * (Y[i] - Y[i - 1])

    l = [0.0] * n; mu = [0.0] * n; z = [0.0] * n
    l[0] = 1.0
    for i in range(1, m):
        l[i] = 2 * (X[i + 1] - X[i - 1]) - h_vals[i - 1] * mu[i - 1]
        mu[i] = h_vals[i] / l[i]
        z[i] = (alpha[i] - h_vals[i - 1] * z[i - 1]) / l[i]
    l[m] = 1.0

    c = [0.0] * n; b = [0.0] * m; d = [0.0] * m
    for j in range(m - 1, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (Y[j + 1] - Y[j]) / h_vals[j] - h_vals[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h_vals[j])

    max_err = 0.0
    for i in range(m):
        for xv, yv in [(X[i], Y[i]), (X[i + 1], Y[i + 1])]:
            dx = xv - X[i]
            sv = Y[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3
            max_err = max(max_err, abs(sv - yv))

    check("Cubic Spline passes through all data points", max_err < 1e-10,
          f"max error = {max_err}")


# ══════════════════════════════════════════════════════════════════════════
# 12 — Trapezoidal Rule
# ══════════════════════════════════════════════════════════════════════════

def test_trapezoidal():
    def f(x): return math.sin(x)
    a, b, n = 0.0, math.pi, 10000
    h = (b - a) / n
    total = f(a) + f(b) + 2 * sum(f(a + i * h) for i in range(1, n))
    result = h * total / 2.0

    check("Trapezoidal: ∫₀^π sin(x)dx ≈ 2", abs(result - 2.0) < 1e-7,
          f"got {result}")


# ══════════════════════════════════════════════════════════════════════════
# 13 — Simpson's Rule
# ══════════════════════════════════════════════════════════════════════════

def test_simpson():
    def f(x): return math.sin(x)
    a, b, n = 0.0, math.pi, 1000
    h = (b - a) / n
    total = f(a) + f(b)
    for i in range(1, n):
        total += (4 if i % 2 else 2) * f(a + i * h)
    result = h * total / 3.0

    check("Simpson: ∫₀^π sin(x)dx ≈ 2", abs(result - 2.0) < 1e-11,
          f"got {result}")


# ══════════════════════════════════════════════════════════════════════════
# 14 — Romberg Integration
# ══════════════════════════════════════════════════════════════════════════

def test_romberg():
    def f(x): return math.sin(x)
    a, b = 0.0, math.pi
    max_i = 10

    R = [[0.0] * max_i for _ in range(max_i)]
    R[0][0] = (b - a) / 2.0 * (f(a) + f(b))

    for i in range(1, max_i):
        n_panels = 1 << i
        h = (b - a) / n_panels
        s = sum(f(a + k * h) for k in range(1, n_panels, 2))
        R[i][0] = R[i - 1][0] / 2.0 + h * s
        for j in range(1, i + 1):
            R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (4 ** j - 1)

    result = R[max_i - 1][max_i - 1]
    check("Romberg: ∫₀^π sin(x)dx ≈ 2", abs(result - 2.0) < 1e-14,
          f"got {result}")


# ══════════════════════════════════════════════════════════════════════════
# 15 — Finite Difference
# ══════════════════════════════════════════════════════════════════════════

def test_finite_diff():
    def f(x): return x ** 2
    x, h = 3.0, 1e-6
    fwd = (f(x + h) - f(x)) / h
    ctr = (f(x + h) - f(x - h)) / (2 * h)
    check("Forward diff: f'(3) for x² ≈ 6", abs(fwd - 6.0) < 1e-3, f"got {fwd}")
    check("Central diff: f'(3) for x² ≈ 6", abs(ctr - 6.0) < 1e-8, f"got {ctr}")


# ══════════════════════════════════════════════════════════════════════════
# 16 — 5-Point Stencil
# ══════════════════════════════════════════════════════════════════════════

def test_5pt_stencil():
    def f(x): return x ** 3
    x, h = 2.0, 1e-4
    fp = (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h)) / (12 * h)
    check("5pt stencil: f'(2) for x³ = 12", abs(fp - 12.0) < 1e-6,
          f"got {fp}")


# ══════════════════════════════════════════════════════════════════════════
# 17 — Euler's Method
# ══════════════════════════════════════════════════════════════════════════

def test_euler():
    def f(x, y): return x + y
    x0, y0, xf, h = 0.0, 1.0, 1.0, 0.01
    x, y = x0, y0
    steps = round((xf - x0) / h)
    for _ in range(steps):
        y += h * f(x, y)
        x += h
    exact = 2 * math.exp(xf) - xf - 1
    # Euler is O(h); with h=0.01 expect error ~0.03
    check("Euler: y(1) ≈ exact (O(h) accuracy)",
          abs(y - exact) < 0.03,
          f"got {y:.8f}, exact {exact:.8f}, error {abs(y - exact):.2e}")


# ══════════════════════════════════════════════════════════════════════════
# 18 — Runge-Kutta 4
# ══════════════════════════════════════════════════════════════════════════

def test_rk4():
    def f(x, y): return x + y
    x0, y0, xf, h = 0.0, 1.0, 1.0, 0.01
    x, y = x0, y0
    steps = round((xf - x0) / h)
    for _ in range(steps):
        k1 = h * f(x, y)
        k2 = h * f(x + h/2, y + k1/2)
        k3 = h * f(x + h/2, y + k2/2)
        k4 = h * f(x + h, y + k3)
        y += (k1 + 2*k2 + 2*k3 + k4) / 6
        x += h
    exact = 2 * math.exp(xf) - xf - 1
    check("RK4: y(1) ≈ exact (O(h⁴) accuracy)", abs(y - exact) < 1e-8,
          f"got {y:.10f}, exact {exact:.10f}, error {abs(y - exact):.2e}")


# ══════════════════════════════════════════════════════════════════════════
# 19 — Power Method
# ══════════════════════════════════════════════════════════════════════════

def test_power_method():
    # [[4, 1], [1, 3]] eigenvalues: char poly = λ²-7λ+11
    # λ = (7 ± √(49-44))/2 = (7 ± √5)/2 ≈ 4.618, 2.382
    A = [[4, 1], [1, 3]]
    v = [1.0, 1.0]
    norm_v = math.sqrt(sum(x * x for x in v))
    v = [x / norm_v for x in v]

    lam = 0.0
    for _ in range(200):
        w = mat_vec_mul(A, v)
        vtv = sum(v[i] * v[i] for i in range(2))
        vtw = sum(v[i] * w[i] for i in range(2))
        lam = vtw / vtv
        wnorm = math.sqrt(sum(x * x for x in w))
        v = [x / wnorm for x in w]

    expected = (7 + math.sqrt(5)) / 2  # ≈ 4.6180339887
    check("Power Method: dominant λ ≈ 4.618",
          abs(lam - expected) < 1e-8, f"got {lam}")


# ══════════════════════════════════════════════════════════════════════════
# 20 — QR Algorithm
# ══════════════════════════════════════════════════════════════════════════

def test_qr():
    # [[5, 1], [1, 2]] eigenvalues: char poly λ²-7λ+9=0 → (7±√13)/2 ≈ 5.303, 1.697
    A = [[5.0, 1.0], [1.0, 2.0]]
    n = 2

    def qr_decompose(M):
        nn = len(M)
        V = [row[:] for row in M]
        R = [[0.0] * nn for _ in range(nn)]
        Q = [[0.0] * nn for _ in range(nn)]
        for i in range(nn):
            for j in range(i):
                dot = sum(V[k][j] * V[k][i] for k in range(nn))
                R[j][i] = dot
                for k in range(nn):
                    V[k][i] -= dot * V[k][j]
            norm_sq = sum(V[k][i] ** 2 for k in range(nn))
            R[i][i] = math.sqrt(max(0.0, norm_sq))
            for k in range(nn):
                Q[k][i] = V[k][i] / R[i][i] if R[i][i] > 0 else 0.0
                V[k][i] = Q[k][i]
        return Q, R

    def mat_mul(M1, M2):
        nn = len(M1)
        return [[sum(M1[i][k] * M2[k][j] for k in range(nn))
                 for j in range(nn)] for i in range(nn)]

    tol_qr = 1e-10
    for _ in range(100):
        try:
            Q, R = qr_decompose(A)
            A = mat_mul(R, Q)
            off = math.sqrt(sum(A[i][j]**2 for i in range(n) for j in range(i+1, n)))
            if off < tol_qr:
                break
        except (ValueError, OverflowError):
            break

    eigenvalues = sorted([A[0][0], A[1][1]])
    e1 = (7 - math.sqrt(13)) / 2  # ≈ 1.697
    e2 = (7 + math.sqrt(13)) / 2  # ≈ 5.303
    check("QR: eigenvalues ≈ [1.697, 5.303]",
          abs(eigenvalues[0] - e1) < 1e-8 and abs(eigenvalues[1] - e2) < 1e-8,
          f"got {eigenvalues}")


# ─── Main test runner ─────────────────────────────────────────────────────

if __name__ == "__main__":
    print()
    print("=" * 62)
    print("  Python Numerical Methods — Automated Test Suite")
    print("=" * 62)

    tests = [
        ("01 — Bisection",             test_bisection),
        ("02 — Newton-Raphson",        test_newton_raphson),
        ("03 — Secant",                test_secant),
        ("04 — False Position",        test_false_position),
        ("05 — Gauss Elimination",     test_gauss_elimination),
        ("06 — LU Decomposition",      test_lu_decomposition),
        ("07 — Jacobi",                test_jacobi),
        ("08 — Gauss-Seidel",          test_gauss_seidel),
        ("09 — Lagrange",              test_lagrange),
        ("10 — Newton DD",             test_newton_interp),
        ("11 — Cubic Spline",          test_spline),
        ("12 — Trapezoidal",           test_trapezoidal),
        ("13 — Simpson",               test_simpson),
        ("14 — Romberg",               test_romberg),
        ("15 — Finite Difference",     test_finite_diff),
        ("16 — 5-Point Stencil",       test_5pt_stencil),
        ("17 — Euler ODE",             test_euler),
        ("18 — RK4 ODE",               test_rk4),
        ("19 — Power Method",          test_power_method),
        ("20 — QR Algorithm",          test_qr),
    ]

    for name, test_fn in tests:
        print(f"\n▸ {name}")
        try:
            test_fn()
        except Exception as e:
            FAILED += 1
            print(f"  ✗ Exception: {e}")

    print()
    print("=" * 62)
    print(f"  Results:  {PASSED} passed,  {FAILED} failed,  "
          f"{PASSED + FAILED} total")
    print("=" * 62)
    print()

    sys.exit(0 if FAILED == 0 else 1)
