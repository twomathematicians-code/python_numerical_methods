"""
utils.py — Shared utilities for the Python Numerical Methods Toolkit.

Provides pretty-print helpers, matrix/vector operations, and mathematical
constants used across all 20 numerical method programs.
"""

import math
from typing import List, Callable, Tuple

# ─── Pretty Print Helpers ──────────────────────────────────────────────────

def print_header(title: str, width: int = 62) -> None:
    """Print a framed header banner."""
    inner = width - 4
    print(f"\n{'╔' + '═' * width + '╗'}")
    print(f"║  {title:<{inner}}║")
    print(f"{'╚' + '═' * width + '╝'}\n")


def print_section(title: str) -> None:
    """Print a section divider."""
    print(f"─── {title} ───")


def print_result(label: str, value: float) -> None:
    """Print a labelled result value."""
    print(f"{label:<32}: {value}")


def print_separator(char: str = "─", width: int = 62) -> None:
    """Print a horizontal separator line."""
    print(char * width)


def print_table_row(*cols: float, widths: List[int] = None) -> None:
    """Print a formatted table row with fixed-width columns."""
    if widths is None:
        widths = [12] * len(cols)
    parts = []
    for val, w in zip(cols, widths):
        parts.append(f"{val:>{w}.8f}")
    print("  ".join(parts))


def print_table_header(*headers: str, widths: List[int] = None) -> None:
    """Print a table header row."""
    if widths is None:
        widths = [12] * len(headers)
    parts = []
    for h, w in zip(headers, widths):
        parts.append(f"{h:>{w}}")
    print("  ".join(parts))


# ─── Matrix / Vector Utilities ─────────────────────────────────────────────

def mat_vec_mul(A: List[List[float]], x: List[float]) -> List[float]:
    """Multiply matrix A by vector x."""
    n = len(A)
    b = [0.0] * n
    for i in range(n):
        for j in range(n):
            b[i] += A[i][j] * x[j]
    return b


def vec_norm(v: List[float]) -> float:
    """Euclidean norm of a vector."""
    return math.sqrt(sum(x * x for x in v))


def vec_max_diff(a: List[float], b: List[float]) -> float:
    """Maximum absolute difference between two vectors."""
    return max(abs(ai - bi) for ai, bi in zip(a, b))


def mat_print(A: List[List[float]], label: str = "Matrix") -> None:
    """Pretty-print a matrix."""
    n = len(A)
    print(f"{label}:")
    for row in A:
        print("  [" + "  ".join(f"{x:12.6f}" for x in row) + " ]")


def vec_print(v: List[float], label: str = "Vector") -> None:
    """Pretty-print a vector."""
    print(f"{label}: [{', '.join(f'{x:.6f}' for x in v)}]")


# ─── Mathematical Helpers ──────────────────────────────────────────────────

def factorial(n: int) -> float:
    """Factorial of n."""
    result = 1.0
    for i in range(2, n + 1):
        result *= i
    return result


def comb(n: int, k: int) -> float:
    """Binomial coefficient C(n, k)."""
    return factorial(n) / (factorial(k) * factorial(n - k))


# ─── Input Helpers ─────────────────────────────────────────────────────────

def prompt_float(prompt: str) -> float:
    """Prompt the user for a float input."""
    return float(input(prompt))


def prompt_int(prompt: str) -> int:
    """Prompt the user for an integer input."""
    return int(input(prompt))
