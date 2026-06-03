#!/usr/bin/env python3
"""
16_higher_order_difference.py — 5-Point Stencil Finite Difference

Implements 5-point stencil formulas for first and second derivatives,
providing O(h^4) accuracy.

f'(x)  ~ [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)] / (12h)         O(h^4)
f''(x) ~ [-f(x+2h) + 16f(x+h) - 30f(x) + 16f(x-h) - f(x-2h)] / (12h^2)  O(h^4)

Run:      python src/16_higher_order_difference.py
"""

import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def f(x: float) -> float:
    return math.exp(x) * math.cos(x)


def fp_exact(x: float) -> float:
    return math.exp(x) * (math.cos(x) - math.sin(x))


def fpp_exact(x: float) -> float:
    return -2.0 * math.exp(x) * math.sin(x)


def main():
    print_header("Higher-Order Central Difference — 5-Point Stencil")

    x = prompt_float("Enter point of evaluation (x): ")
    n_steps = prompt_int("Enter number of step-size refinements: ")

    print(f"\nExact f'({x})   = {fp_exact(x):.10f}")
    print(f"Exact f''({x})  = {fpp_exact(x):.10f}\n")

    header = (f"{'Step':>4}  {'h':>10}  {'5pt fp':>14}  {'Err_fp':>10}"
              f"  {'5pt fpp':>14}  {'Err_fpp':>10}")
    print(header)
    print_separator(width=72)

    h = 0.5
    for i in range(n_steps):
        # 5-point first derivative
        fp = (-f(x + 2*h) + 8.0*f(x + h) - 8.0*f(x - h) + f(x - 2*h)) / (12.0 * h)
        err_fp = abs(fp - fp_exact(x))

        # 5-point second derivative
        fpp = (-f(x + 2*h) + 16.0*f(x + h) - 30.0*f(x)
              + 16.0*f(x - h) - f(x - 2*h)) / (12.0 * h * h)
        err_fpp = abs(fpp - fpp_exact(x))

        print(f"{i:>4}  {h:>10.2e}  {fp:>14.10f}  {err_fp:>10.4e}"
              f"  {fpp:>14.10f}  {err_fpp:>10.4e}")
        h /= 2.0

    print_separator(width=72)
    return 0


if __name__ == "__main__":
    main()
