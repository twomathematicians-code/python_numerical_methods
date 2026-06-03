#!/usr/bin/env python3
"""
05_gauss_elimination.py — Gaussian Elimination with Partial Pivoting

Solves Ax = b by transforming the augmented matrix [A|b] into upper
triangular form via row operations, then back-substituting.

Time complexity: O(n^3)

Run:      python src/05_gauss_elimination.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, prompt_int


def main():
    print_header("Gaussian Elimination — Linear System Solver")

    n = prompt_int("Enter the number of equations (n): ")

    print(f"\nEnter the augmented matrix [A|b] row by row:")
    A = []
    for i in range(n):
        row = []
        for j in range(n + 1):
            row.append(float(input(f"  A[{i}][{j}]: ")))
        A.append(row)

    # ── Forward elimination with partial pivoting ──
    for col in range(n):
        pivot = col
        for row in range(col + 1, n):
            if abs(A[row][col]) > abs(A[pivot][col]):
                pivot = row

        if abs(A[pivot][col]) < 1e-12:
            print("Error: Matrix is singular or nearly singular.")
            return 1

        A[col], A[pivot] = A[pivot], A[col]

        for row in range(col + 1, n):
            factor = A[row][col] / A[col][col]
            for j in range(col, n + 1):
                A[row][j] -= factor * A[col][j]

    # ── Back substitution ──
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = A[i][n]
        for j in range(i + 1, n):
            x[i] -= A[i][j] * x[j]
        x[i] /= A[i][i]

    print("\nSolution vector:")
    for i in range(n):
        print(f"  x[{i}] = {x[i]:.10f}")

    return 0


if __name__ == "__main__":
    main()
