#!/usr/bin/env python3
"""
10_newton_interpolation.py — Newton's Divided Difference Interpolation

Constructs the interpolating polynomial using divided differences,
allowing incremental addition of new data points.

P(x) = f[x0] + f[x0,x1](x-x0) + f[x0,x1,x2](x-x0)(x-x1) + ...

Run:      python src/10_newton_interpolation.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def main():
    print_header("Newton's Divided Difference Interpolation")

    n = prompt_int("Enter number of data points: ")

    X, Y = [], []
    print("\nEnter data points (x, y):")
    for i in range(n):
        X.append(float(input(f"  Point {i} — x: ")))
        Y.append(float(input(f"             y: ")))

    xval = prompt_float("\nEnter value to interpolate at (x): ")

    # ── Build divided-difference table ──
    dd = [[0.0] * n for _ in range(n)]
    for i in range(n):
        dd[i][0] = Y[i]

    for j in range(1, n):
        for i in range(n - j):
            dd[i][j] = (dd[i + 1][j - 1] - dd[i][j - 1]) / (X[i + j] - X[i])

    # Display the table
    print("\nDivided Difference Table:")
    for i in range(n):
        row_str = ""
        for j in range(n - i):
            row_str += f"{dd[i][j]:>14.6f}"
        print(row_str)

    # ── Evaluate Newton polynomial ──
    result = dd[0][0]
    product = 1.0

    poly_terms = f"{dd[0][0]:.8f}"
    for j in range(1, n):
        product *= (xval - X[j - 1])
        result += dd[0][j] * product
        poly_terms += f" + ({dd[0][j]:.8f})"
        for k in range(j):
            poly_terms += f" * (x - {X[k]})"

    print(f"\nPolynomial terms:\n  {poly_terms}")
    print_separator()
    print(f"P({xval}) = {result:.10f}")

    return 0


if __name__ == "__main__":
    main()
