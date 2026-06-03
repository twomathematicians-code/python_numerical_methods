#!/usr/bin/env python3
"""
06_lu_decomposition.py — LU Decomposition (Doolittle's Method)

Factors A = L * U, where L is lower-triangular with unit diagonal and U is
upper-triangular.  Then solves Ax = b as:
  1) Forward substitution:  Ly = b
  2) Back substitution:      Ux = y

Run:      python src/06_lu_decomposition.py
"""

import sys, os, copy
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import (print_header, prompt_int, mat_print, vec_print,
                   mat_vec_mul)


def main():
    print_header("LU Decomposition (Doolittle) — Linear System Solver")

    n = prompt_int("Enter the matrix dimension (n): ")

    print("\nEnter matrix A row by row:")
    A = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(float(input(f"  A[{i}][{j}]: ")))
        A.append(row)

    print("\nEnter vector b:")
    b = []
    for i in range(n):
        b.append(float(input(f"  b[{i}]: ")))

    # ── LU factorisation (Doolittle) ──
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

    mat_print(L, "\nL matrix")
    mat_print(U, "\nU matrix")

    # ── Forward substitution: Ly = b ──
    y = [0.0] * n
    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))

    # ── Back substitution: Ux = y ──
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]

    print("\nSolution vector:")
    for i in range(n):
        print(f"  x[{i}] = {x[i]:.10f}")

    # Verification
    Ax = mat_vec_mul(A, x)
    print("\nVerification (Ax vs b):")
    for i in range(n):
        print(f"  Ax[{i}] = {Ax[i]:.10f}  |  b[{i}] = {b[i]:.10f}"
              f"  |  error = {abs(Ax[i] - b[i]):.2e}")

    return 0


if __name__ == "__main__":
    main()
