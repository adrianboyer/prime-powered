"""Irreducible polynomials in unicritically generated semigroups.

Computational side of Section 3 of arXiv:2510.10310.  The objects:
TBA
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterator

from sympy import Poly, primefactors, symbols

from .dynamics import prime_power_forms

x = symbols("x")

SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71)


def unicritical(c: int, d: int) -> Poly:
    return Poly(x**d + c, x)
