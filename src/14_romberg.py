#!/usr/bin/env python3
"""
14_romberg.py — Romberg Integration

Applies Richardson extrapolation to the trapezoidal rule, producing
increasingly accurate estimates.

R(0,0) = Trapezoidal with 1 panel
R(i,0) = Trapezoidal with 2^i panels
R(i,j) = R(i,j-1) + [R(i,j-1) - R(i-1,j-1)] / (4^j - 1)

Run:      python src/14_romberg.py
"""

import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def f(x: float) -> float:
    return math.sin(x)


def main():
    print_header("Romberg Integration — Richardson Extrapolation")

    a = prompt_float("Enter lower limit (a): ")
    b = prompt_float("Enter upper limit (b): ")
    max_iter = prompt_int("Enter maximum Romberg iterations (typically 5-10): ")
    tol = prompt_float("Enter tolerance: ")

    R = [[0.0] * max_iter for _ in range(max_iter)]
    R[0][0] = (b - a) / 2.0 * (f(a) + f(b))

    print("\nRomberg Extrapolation Table:\n")

    for i in range(1, max_iter):
        n_panels = 1 << i  # 2^i
        h = (b - a) / n_panels

        # Composite trapezoidal with 2^i panels (only new odd points)
        s = sum(f(a + k * h) for k in range(1, n_panels, 2))
        R[i][0] = R[i - 1][0] / 2.0 + h * s

        # Richardson extrapolation
        for j in range(1, i + 1):
            factor = 4.0 ** j
            R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (factor - 1.0)

        row_str = f"n={n_panels:>6}"
        for j in range(i + 1):
            row_str += f"  {R[i][j]:>18.10f}"
        print(row_str)

        if abs(R[i][i] - R[i - 1][i - 1]) < tol:
            print(f"\nConverged at iteration {i}.")
            break

    best = R[max_iter - 1][max_iter - 1]
    print(f"\nBest estimate = {best:.10f}")

    return 0


if __name__ == "__main__":
    main()
