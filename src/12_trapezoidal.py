#!/usr/bin/env python3
"""
12_trapezoidal.py — Composite Trapezoidal Rule for Numerical Integration

Approximates integral f(x) dx from a to b using n subintervals.

integral = (h/2) * [f(x0) + 2*sum(f(x_i)) + f(x_n)]
where h = (b - a) / n

Run:      python src/12_trapezoidal.py
"""

import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def f(x: float) -> float:
    # Example: f(x) = sin(x), integral from 0 to pi = 2
    return math.sin(x)


def main():
    print_header("Composite Trapezoidal Rule — Numerical Integration")

    a = prompt_float("Enter lower limit (a): ")
    b = prompt_float("Enter upper limit (b): ")
    n = prompt_int("Enter number of subintervals (n): ")

    if n <= 0:
        print("Error: n must be positive.")
        return 1

    h = (b - a) / n

    print()
    header = f"{'i':>4}  {'x_i':>16}  {'f(x_i)':>16}  {'weight':>8}"
    print(header)
    print_separator(width=56)

    total = 0.0
    for i in range(n + 1):
        xi = a + i * h
        fi = f(xi)
        weight = 1.0 if (i == 0 or i == n) else 2.0
        total += weight * fi
        print(f"{i:>4}  {xi:>16.10f}  {fi:>16.10f}  {weight:>8.0f}")

    result = h * total / 2.0

    print_separator(width=56)
    print(f"\nh = {h}")
    print(f"Approximate integral = {result:.10f}")
    print(f"(Error bound approx {-(b - a) * h ** 2 / 12.0:.4e})")

    return 0


if __name__ == "__main__":
    main()
