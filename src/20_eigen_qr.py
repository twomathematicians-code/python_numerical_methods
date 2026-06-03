#!/usr/bin/env python3
"""
20_eigen_qr.py — QR Algorithm for All Eigenvalues

Iteratively decomposes A = QR, then replaces A with RQ.  Under mild
conditions this converges to an upper-triangular matrix whose diagonal
contains the eigenvalues.

Uses modified Gram-Schmidt for QR factorisation (more numerically stable).

Run:      python src/20_eigen_qr.py
"""

import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_int, prompt_float


def qr_decompose(A):
    """Modified Gram-Schmidt QR factorisation (column-wise). Returns Q, R."""
    n = len(A)
    V = [row[:] for row in A]
    R = [[0.0] * n for _ in range(n)]
    Q = [[0.0] * n for _ in range(n)]

    for i in range(n):
        # Orthogonalise column i against all previous columns j
        for j in range(i):
            dot = sum(V[k][j] * V[k][i] for k in range(n))
            R[j][i] = dot
            for k in range(n):
                V[k][i] -= dot * V[k][j]

        # Column norm
        norm_sq = sum(V[k][i] ** 2 for k in range(n))
        R[i][i] = math.sqrt(max(0.0, norm_sq))

        # Normalise into Q, and update V to store normalised columns
        for k in range(n):
            Q[k][i] = V[k][i] / R[i][i] if R[i][i] > 0 else 0.0
            V[k][i] = Q[k][i]

    return Q, R


def mat_mul(A, B):
    """Matrix multiplication C = A * B."""
    n = len(A)
    C = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C


def main():
    print_header("QR Algorithm — All Eigenvalues")

    n = prompt_int("Enter matrix dimension (n): ")

    print("\nEnter matrix A row by row:")
    A = []
    for i in range(n):
        row = [float(input(f"  A[{i}][{j}]: ")) for j in range(n)]
        A.append(row)

    max_iter = prompt_int("Enter maximum QR iterations: ")
    tol = prompt_float("Enter tolerance: ")

    print("\nQR iteration progress:")
    print_separator()

    for iteration in range(max_iter):
        Q, R = qr_decompose(A)
        A = mat_mul(R, Q)  # A_{k+1} = R_k * Q_k

        # Check off-diagonal convergence
        off = math.sqrt(sum(A[i][j] ** 2 for i in range(n) for j in range(i + 1, n)))

        if (iteration + 1) % 10 == 0 or iteration == 0 or off < tol:
            diag = ", ".join(f"{A[i][i]:.6f}" for i in range(n))
            print(f"Iter {iteration + 1:>4}  off-diag = {off:.6e}  diagonal: [{diag}]")

        if off < tol:
            print_separator()
            print(f"\nConverged after {iteration + 1} iterations.")
            break

    print("\nApproximate eigenvalues (diagonal of converged matrix):")
    for i in range(n):
        print(f"  lambda_{i} = {A[i][i]:.10f}")

    return 0


if __name__ == "__main__":
    main()
