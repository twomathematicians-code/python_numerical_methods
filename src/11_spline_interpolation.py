#!/usr/bin/env python3
"""
11_spline_interpolation.py — Natural Cubic Spline Interpolation

Fits a smooth piecewise cubic polynomial through n+1 data points.
Natural boundary condition: second derivative = 0 at both endpoints.

S_i(x) = a + b*(x-xi) + c*(x-xi)^2 + d*(x-xi)^3

Run:      python src/11_spline_interpolation.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def main():
    print_header("Natural Cubic Spline Interpolation")

    n = prompt_int("Enter number of data points: ")

    X, Y = [], []
    print("\nEnter data points (x, y), sorted by x:")
    for i in range(n):
        X.append(float(input(f"  Point {i} — x: ")))
        Y.append(float(input(f"             y: ")))

    xval = prompt_float("\nEnter value to interpolate at (x): ")

    m = n - 1
    h = [X[i + 1] - X[i] for i in range(m)]

    # ── Solve tridiagonal system for second derivatives ──
    alpha = [0.0] * m
    for i in range(1, m):
        alpha[i] = (3.0 / h[i]) * (Y[i + 1] - Y[i]) - \
                   (3.0 / h[i - 1]) * (Y[i] - Y[i - 1])

    l = [0.0] * n
    mu = [0.0] * n
    z = [0.0] * n
    l[0] = 1.0

    for i in range(1, m):
        l[i] = 2.0 * (X[i + 1] - X[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    l[m] = 1.0

    c = [0.0] * n
    b = [0.0] * m
    d = [0.0] * m

    for j in range(m - 1, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (Y[j + 1] - Y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j])

    # Display spline coefficients
    print("\nSpline coefficients (S_i(x) = a + b(x-xi) + c(x-xi)^2 + d(x-xi)^3):")
    for i in range(m):
        print(f"  S{i} on [{X[i]}, {X[i + 1]}]:")
        print(f"    a={Y[i]:.6f}  b={b[i]:.6f}  c={c[i]:.6f}  d={d[i]:.6f}")

    # ── Evaluate spline at xval ──
    idx = 0
    for i in range(m):
        if X[i] <= xval <= X[i + 1]:
            idx = i
            break
    else:
        idx = m - 1

    dx = xval - X[idx]
    result = Y[idx] + b[idx] * dx + c[idx] * dx ** 2 + d[idx] * dx ** 3

    print_separator()
    print(f"S({xval}) = {result:.10f}")

    return 0


if __name__ == "__main__":
    main()
