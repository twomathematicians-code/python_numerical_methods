#!/usr/bin/env python3
"""
08_gauss_seidel.py — Gauss-Seidel Iterative Method for Linear Systems

Similar to Jacobi, but uses the latest updated values immediately (in-place),
which typically leads to faster convergence.

x_i^{(k+1)} = (1 / a_{ii}) * (b_i - sum_{j<i} a_{ij}*x_j^{(k+1)} - sum_{j>i} a_{ij}*x_j^{(k)})

Run:      python src/08_gauss_seidel.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_int, prompt_float


def main():
    print_header("Gauss-Seidel Iterative Method — Linear System Solver")

    n = prompt_int("Enter the matrix dimension (n): ")

    print("\nEnter matrix A row by row:")
    A = []
    for i in range(n):
        row = [float(input(f"  A[{i}][{j}]: ")) for j in range(n)]
        A.append(row)

    print("\nEnter vector b:")
    b = [float(input(f"  b[{i}]: ")) for i in range(n)]

    tol = prompt_float("Enter tolerance: ")
    max_iter = prompt_int("Enter maximum iterations: ")

    x = [0.0] * n

    print()
    print(f"{'Iter':>6}  {'Error':>14}  x vector")
    print_separator()

    for iteration in range(max_iter):
        x_old = x[:]
        for i in range(n):
            s = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (b[i] - s) / A[i][i]

        err = max(abs(x[i] - x_old[i]) for i in range(n))

        vec_str = "[" + ", ".join(f"{v:.6f}" for v in x) + "]"
        print(f"{iteration + 1:>6}  {err:>14.4e}  {vec_str}")

        if err < tol:
            print_separator()
            print(f"\nConverged after {iteration + 1} iterations.")
            break

    print("\nFinal solution:")
    for i in range(n):
        print(f"  x[{i}] = {x[i]:.10f}")

    return 0


if __name__ == "__main__":
    main()
