#!/usr/bin/env python3
"""
13_simpson.py — Composite Simpson's 1/3 Rule for Numerical Integration

Uses quadratic polynomial approximations over pairs of subintervals.
Requires an even number of subintervals.

integral = (h/3) * [f(x0) + 4*sum(f(x_odd)) + 2*sum(f(x_even)) + f(x_n)]

Run:      python src/13_simpson.py
"""

import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def f(x: float) -> float:
    return math.sin(x)


def main():
    print_header("Composite Simpson's 1/3 Rule — Numerical Integration")

    a = prompt_float("Enter lower limit (a): ")
    b = prompt_float("Enter upper limit (b): ")
    n = prompt_int("Enter number of subintervals (n, must be even): ")

    if n <= 0 or n % 2 != 0:
        print("Error: n must be a positive even number.")
        return 1

    h = (b - a) / n

    total = f(a) + f(b)

    header = f"{'i':>4}  {'x_i':>16}  {'f(x_i)':>16}  {'weight':>8}"
    print()
    print(header)
    print_separator(width=56)

    print(f"{0:>4}  {a:>16.10f}  {f(a):>16.10f}  {1:>8}")

    for i in range(1, n):
        xi = a + i * h
        fi = f(xi)
        weight = 4.0 if i % 2 != 0 else 2.0
        total += weight * fi
        print(f"{i:>4}  {xi:>16.10f}  {fi:>16.10f}  {weight:>8.0f}")

    print(f"{n:>4}  {b:>16.10f}  {f(b):>16.10f}  {1:>8}")

    result = h * total / 3.0

    print_separator(width=56)
    print(f"\nh = {h}")
    print(f"Approximate integral = {result:.10f}")
    print(f"(Error bound approx {-(b - a) * h ** 4 / 180.0:.4e})")

    return 0


if __name__ == "__main__":
    main()
