"""Irreducible polynomials in unicritically generated semigroups.

Computational side of Section 3 of arXiv:2510.10310.  The objects:

* f_i(x) = x^d + c_i, and the semigroup G = ⟨f_1, …, f_s⟩ under composition.
  A *word* (i_1, …, i_n) stands for f_{i_1} ∘ f_{i_2} ∘ ⋯ ∘ f_{i_n}; G is free
  (Prop. 3.9), so distinct words are distinct polynomials and "length" is
  well defined.

* Prop. 3.1 — if w is irreducible (of even degree when d is even) and
  w ∘ (x^d + c) is reducible, then w(c) is a p-th power for some prime p | d.
  ``prop_3_1_witnesses`` searches for reducible compositions and checks the
  necessary condition on every one.

* Prop. 3.2 — x^d + c irreducible over Q  ⟹  every iterate is irreducible
  (stability).  ``is_stable`` checks it to a given depth.

* Props. 3.5 / 3.7 / 3.8 — explicit infinite families of irreducible
  polynomials f_1^3 ∘ g, f_1^3 ∘ f_2 ∘ f_1 ∘ g, f_1^3 ∘ f_2^3 ∘ g (f_1^4 ∘ g when
  d = 2), whose existence gives the positive proportion in Theorem 1.1.
  ``family_check`` verifies a family on all g up to a given length.

Irreducibility over Q is decided two ways.  ``is_irreducible`` is SymPy's
exact test.  ``certificate`` looks for a prime p such that the polynomial is
irreducible mod p — a one-line *proof* of irreducibility over Q that costs
almost nothing even when the coefficients have hundreds of digits (as they
do in the exceptional semigroups of Theorem 1.1).  A missing certificate
proves nothing: x⁴ + 1 is irreducible over Q but reducible mod every prime.
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


def compose(outer: Poly, inner: Poly) -> Poly:
    """outer ∘ inner."""
    return outer.compose(inner)


def word_polynomial(word: tuple[int, ...], coeffs, d: int) -> Poly:
    """f_{w_1} ∘ f_{w_2} ∘ ⋯ ∘ f_{w_n}  for f_i = x^d + coeffs[i]  (the empty word is x)."""
    poly = Poly(x, x)
    for i in reversed(word):
        poly = poly**d + coeffs[i]
    return poly


def words(s: int, max_len: int, min_len: int = 1) -> Iterator[tuple[int, ...]]:
    for n in range(min_len, max_len + 1):
        yield from product(range(s), repeat=n)


# --------------------------------------------------------------------------
# irreducibility
# --------------------------------------------------------------------------

def is_irreducible(poly: Poly) -> bool:
    """Exact irreducibility over Q (Gauss: for a monic integer polynomial, same as over Z)."""
    return bool(poly.is_irreducible)


def irreducible_mod(poly: Poly, p: int) -> bool:
    """Irreducible over F_p?  Only meaningful when p does not divide the leading coefficient."""
    return bool(Poly(poly.as_expr(), x, modulus=p).is_irreducible)


def certificate(poly: Poly, primes=SMALL_PRIMES) -> int | None:
    """A prime p with poly irreducible mod p, if one of ``primes`` works.  Such a p proves
    irreducibility over Q (for monic polynomials, which all of ours are)."""
    lead = int(poly.LC())
    for p in primes:
        if lead % p == 0:
            continue
        if irreducible_mod(poly, p):
            return p
    return None


def decide(poly: Poly, primes=SMALL_PRIMES) -> tuple[bool, str]:
    """(irreducible?, how): cheap certificate first, exact factorisation as the fallback."""
    p = certificate(poly, primes)
    if p is not None:
        return True, f"irreducible mod {p}"
    return is_irreducible(poly), "exact factorisation"


# --------------------------------------------------------------------------
# Propositions 3.1 and 3.2
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class CompositionWitness:
    a: int          # w = x^d + a
    b: int          # u = x^d + b
    reducible: bool
    w_of_u0: int    # w(u(0)) = b^d + a
    pth_power: bool # is w(u(0)) = y^p for some prime p | d ?


def prop_3_1_witnesses(d: int, a_values, b_values) -> list[CompositionWitness]:
    """For w = x^d + a irreducible and u = x^d + b, record whether w ∘ u is reducible and whether
    w(u(0)) is a p-th power.  Prop. 3.1: reducible ⟹ p-th power (never the other way round)."""
    out = []
    for a in a_values:
        w = unicritical(a, d)
        if not is_irreducible(w):
            continue
        for b in b_values:
            u = unicritical(b, d)
            v = b**d + a
            forms = [f for f in prime_power_forms(v, d) if f[0] == 1 or f[2] % 2 == 1]  # y^p, y ∈ Z
            out.append(CompositionWitness(a, b, not is_irreducible(compose(w, u)), v, bool(forms)))
    return out


def is_stable(c: int, d: int, depth: int) -> bool:
    """Prop. 3.2 at finite depth: f, f², …, f^depth all irreducible (f = x^d + c)."""
    poly = unicritical(c, d)
    for _ in range(depth):
        if not decide(poly)[0]:
            return False
        poly = poly**d + c
    return True


# --------------------------------------------------------------------------
# powered fixed points and 2-cycles (Definition 3.4)
# --------------------------------------------------------------------------

def powered_fixed_points(c: int, d: int) -> list[tuple[int, int]]:
    """(y, p) with p | d prime and f(y^p) = y^p, i.e. c = y^p − y^{pd}.  Exhaustive: |y| is tiny."""
    out = []
    for p in primefactors(d):
        y = 0
        while True:
            for s in ((y,) if y == 0 else (y, -y)):
                if s**p - s**(p * d) == c:
                    out.append((s, p))
            if y >= 1 and abs(y**(p * d)) > 2 * abs(c) + 2:
                break
            y += 1
    return out


def powered_two_cycles(c: int, d: int) -> list[tuple[int, int]]:
    """(y, p) with f(f(y^p)) = y^p and f(y^p) ≠ y^p."""
    out = []
    for p in primefactors(d):
        y = 0
        while True:
            for s in ((y,) if y == 0 else (y, -y)):
                v = s**p
                f1 = v**d + c
                if f1 != v and f1**d + c == v:
                    out.append((s, p))
            if y >= 1 and abs(y**(p * d)) > 2 * abs(c) + 2:
                break
            y += 1
    return out


# --------------------------------------------------------------------------
# semigroups
# --------------------------------------------------------------------------

def semigroup_report(coeffs, d: int, max_len: int, exact_up_to_degree: int = 200) -> dict:
    """Irreducibility of every element of ⟨x^d + c_1, …, x^d + c_s⟩ up to a given word length.

    Polynomials of degree ≤ ``exact_up_to_degree`` are decided exactly; larger ones are only
    *certified* irreducible when a mod-p certificate exists, otherwise reported as undecided.
    """
    coeffs = tuple(coeffs)
    s = len(coeffs)
    by_length = {}
    for n in range(1, max_len + 1):
        rows = []
        for w in product(range(s), repeat=n):
            poly = word_polynomial(w, coeffs, d)
            p = certificate(poly)
            if p is not None:
                rows.append((w, True, f"mod {p}"))
            elif poly.degree() <= exact_up_to_degree:
                rows.append((w, is_irreducible(poly), "exact"))
            else:
                rows.append((w, None, "undecided"))
        irreducible = sum(1 for _, r, _ in rows if r is True)
        reducible = sum(1 for _, r, _ in rows if r is False)
        by_length[n] = {
            "words": len(rows),
            "irreducible": irreducible,
            "reducible": reducible,
            "undecided": len(rows) - irreducible - reducible,
            "reducible_words": [w for w, r, _ in rows if r is False],
            "undecided_words": [w for w, r, _ in rows if r is None],
        }
    total = sum(v["words"] for v in by_length.values())
    irr = sum(v["irreducible"] for v in by_length.values())
    return {"coeffs": coeffs, "d": d, "by_length": by_length,
            "proportion_irreducible_lower_bound": irr / total}


def is_free_up_to(coeffs, d: int, max_len: int) -> bool:
    """Prop. 3.9 at finite length: distinct words give distinct polynomials."""
    seen = set()
    for w in words(len(coeffs), max_len):
        key = tuple(int(t) for t in word_polynomial(w, coeffs, d).all_coeffs())
        if key in seen:
            return False
        seen.add(key)
    return True


def family_check(prefix: tuple[int, ...], coeffs, d: int, max_len_g: int) -> dict:
    """Check that F ∘ g is irreducible for every g ∈ G of length ≤ max_len_g (and g = x),
    where F is the word ``prefix``.  E.g. prefix = (0, 0, 0) is f_1^3 (Prop. 3.5)."""
    results = {}
    for g in [()] + list(words(len(coeffs), max_len_g)):
        poly = word_polynomial(prefix + g, coeffs, d)
        results[g] = decide(poly)
    return {"prefix": prefix, "all_irreducible": all(r[0] for r in results.values()), "results": results}
