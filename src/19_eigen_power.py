#!/usr/bin/env python3
"""
19_eigen_power.py — Power Method for Dominant Eigenvalue

Iteratively applies A*v and normalises to converge to the eigenvector
corresponding to the eigenvalue of largest magnitude.

v_{k+1} = A * v_k / ||A * v_k||
lambda  ~ v_k^T * A * v_k / v_k^T * v_k   (Rayleigh quotient)

Run:      python src/19_eigen_power.py
"""

import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import (print_header, print_separator, prompt_int, prompt_float,
                   mat_vec_mul, vec_norm)


def main():
    print_header("Power Method — Dominant Eigenvalue")

    n = prompt_int("Enter matrix dimension (n): ")

    print("\nEnter matrix A row by row:")
    A = []
    for i in range(n):
        row = [float(input(f"  A[{i}][{j}]: ")) for j in range(n)]
        A.append(row)

    tol = prompt_float("Enter tolerance: ")
    max_iter = prompt_int("Enter maximum iterations: ")

    # Initial eigenvector guess (all 1s, normalised)
    v = [1.0] * n
    norm = vec_norm(v)
    v = [x / norm for x in v]

    print()
    print(f"{'Iter':>6}  {'Eigenvalue':>16}  {'Delta':>16}  Eigenvector")
    print_separator()

    lambda_old = 0.0
    for iteration in range(max_iter):
        w = mat_vec_mul(A, v)

        # Rayleigh quotient
        vtv = sum(vi * vi for vi in v)
        vtw = sum(vi * wi for vi, wi in zip(v, w))
        lam = vtw / vtv

        # Normalise w
        wnorm = vec_norm(w)
        w = [x / wnorm for x in w]

        delta = abs(lam - lambda_old)
        vec_str = "[" + ", ".join(f"{x:.6f}" for x in w) + "]"
        print(f"{iteration + 1:>6}  {lam:>16.10f}  {delta:>16.4e}  {vec_str}")

        if delta < tol:
            print_separator()
            print(f"\nConverged after {iteration + 1} iterations.")
            print(f"Dominant eigenvalue lambda = {lam:.10f}")
            print(f"Corresponding eigenvector:\n  {vec_str}")
            return 0

        v = w
        lambda_old = lam

    print(f"\nMethod did not converge within {max_iter} iterations.")
    return 0


if __name__ == "__main__":
    main()
