#!/usr/bin/env python3
"""
09_lagrange_interpolation.py — Lagrange Polynomial Interpolation

Given (n+1) data points, constructs the unique degree-n polynomial that
passes through all points using Lagrange basis polynomials:

P(x) = sum_{i=0}^{n} y_i * L_i(x)
where L_i(x) = product_{j!=i} (x - x_j) / (x_i - x_j)

Run:      python src/09_lagrange_interpolation.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def main():
    print_header("Lagrange Polynomial Interpolation")

    n = prompt_int("Enter number of data points: ")

    X, Y = [], []
    print("\nEnter data points (x, y):")
    for i in range(n):
        X.append(float(input(f"  Point {i} — x: ")))
        Y.append(float(input(f"             y: ")))

    xval = prompt_float("\nEnter value to interpolate at (x): ")

    # ── Compute Lagrange interpolation ──
    result = 0.0
    print(f"\nLagrange basis values at x = {xval}:")
    for i in range(n):
        Li = 1.0
        for j in range(n):
            if j != i:
                Li *= (xval - X[j]) / (X[i] - X[j])
        term = Li * Y[i]
        result += term
        print(f"  L{i}({xval}) = {Li:.10f}  -> term = {term:.10f}")

    print_separator()
    print(f"P({xval}) = {result:.10f}")

    return 0


if __name__ == "__main__":
    main()
