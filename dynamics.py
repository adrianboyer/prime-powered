"""Prime-powered images of the unicritical polynomials f(x) = x^d + c over Z.

This module is the computational side of Theorem 2.1 of

    Bhardwaj, Boyer-Paulet, Hindes, Qiu, Sun,
    "Prime-powered images and irreducible polynomials in dynamical semigroups",
    arXiv:2510.10310 (accepted, Proc. Amer. Math. Soc. Ser. B).

The theorem:  let f(x) = x^d + c with c ∈ Z nonzero and d ≥ 2, and let
N = 4 if d = 2 and N = 3 if d ≥ 3.  If f^N(α) = ε·y^p for some integers
α, y, some sign ε and some prime p | d, then α is preperiodic and ε·y^p is
periodic.  More precisely one of four explicit descriptions holds:

    (1) d = 2:                 α = ±ε y²,  ε y² is fixed or has exact period 2
    (2) d ≥ 3 odd:             α = ε y^p,  ε y^p is a fixed point
    (3) d ≥ 4 even, c ≠ −1:    α = ±ε y^p, ε y^p is a fixed point
    (4) d ≥ 4 even, c = −1:    f(α) = ε y^p, ε y^p ∈ {0, −1} has exact period 2

NOTE ON VERSIONS.  Statement (4) is quoted above in the form it takes in the
ACCEPTED version.  The public arXiv v1 (2510.10310) prints it as α = ±ε y^p,
which is the d = 2 relation and cannot hold for even d ≥ 4; this was caught and
corrected before acceptance.  ``statement_as_printed`` below keeps the arXiv v1
reading so the difference stays visible to anyone reading the preprint.

and N is sharp (Remark 1.5): f(x) = x² − 460 has f³(22) = 114² with 22 *not*
preperiodic, and f(x) = x^d − r^d has f²(r) = −r^d.

Nothing here assumes the theorem.  ``search_theorem`` finds every solution of
f^N(α) = ε y^p in a window by brute force and ``classify`` checks each one
against the four descriptions independently, by computing orbits.  The
sharpness search ``search_near_misses`` looks one iterate earlier and finds
the counterexamples that show N cannot be lowered — including the two the
paper quotes and the ones it doesn't.
"""

from __future__ import annotations

from dataclasses import dataclass

from sympy import integer_nthroot, primefactors


# --------------------------------------------------------------------------
# iteration and orbits
# --------------------------------------------------------------------------

def apply(alpha: int, c: int, d: int) -> int:
    return alpha**d + c


def iterate(alpha: int, c: int, d: int, n: int) -> int:
    """f^n(α) for f(x) = x^d + c, exact integer arithmetic."""
    for _ in range(n):
        alpha = alpha**d + c
    return alpha


def escapes(x: int, c: int, d: int) -> bool:
    """True if the forward orbit of x provably goes to infinity.

    If |x| ≥ 2 and |x|^d − |x| − |c| > 0 then |f(x)| ≥ |x|^d − |c| > |x| ≥ 2, and
    g(t) = t^d − t − |c| is increasing for t ≥ 1, so the same inequality holds
    for f(x): the orbit increases strictly forever.
    """
    ax = abs(x)
    return ax >= 2 and ax**d - ax - abs(c) > 0


@dataclass(frozen=True)
class Orbit:
    """Preperiodicity data of an integer point.  tail = 0 means periodic."""

    preperiodic: bool
    tail: int | None      # number of steps before the orbit enters its cycle
    period: int | None    # exact period of the cycle
    points: tuple[int, ...]


def orbit(alpha: int, c: int, d: int, max_steps: int = 10_000) -> Orbit:
    """Decide whether α is preperiodic under x^d + c and, if so, find its tail and period."""
    seen: dict[int, int] = {}
    points: list[int] = []
    x = alpha
    for step in range(max_steps):
        if x in seen:
            first = seen[x]
            return Orbit(True, first, step - first, tuple(points))
        if escapes(x, c, d):
            return Orbit(False, None, None, tuple(points))
        seen[x] = step
        points.append(x)
        x = x**d + c
    raise RuntimeError("orbit undecided after max_steps — should not happen for integer unicritical maps")


def is_preperiodic(alpha: int, c: int, d: int) -> bool:
    return orbit(alpha, c, d).preperiodic


def is_periodic(alpha: int, c: int, d: int) -> bool:
    o = orbit(alpha, c, d)
    return o.preperiodic and o.tail == 0


def exact_period(alpha: int, c: int, d: int) -> int | None:
    o = orbit(alpha, c, d)
    return o.period if (o.preperiodic and o.tail == 0) else None


# --------------------------------------------------------------------------
# prime powers
# --------------------------------------------------------------------------

def prime_power_forms(v: int, d: int) -> list[tuple[int, int, int]]:
    """All ways to write v = ε·y^p with ε = ±1, y ≥ 0 an integer and p a prime dividing d.

    Returned as (ε, y, p).  v = 0 is 0^p for every p (y = 0 is allowed in the theorem).
    """
    out = []
    for p in primefactors(d):
        if v == 0:
            out.append((1, 0, p))
            continue
        root, exact = integer_nthroot(abs(v), p)
        if exact:
            out.append((1 if v > 0 else -1, int(root), p))
    return out


# --------------------------------------------------------------------------
# Theorem 2.1: search and classify
# --------------------------------------------------------------------------

def threshold_iterate(d: int) -> int:
    """N in Theorem 2.1: 4 when d = 2, 3 when d ≥ 3."""
    return 4 if d == 2 else 3


@dataclass(frozen=True)
class Solution:
    c: int
    d: int
    alpha: int
    n: int
    value: int          # f^n(α)
    eps: int
    y: int
    p: int
    alpha_orbit: Orbit
    value_orbit: Orbit


def statement(sol: Solution) -> tuple[bool, str]:
    """Check a solution of f^N(α) = ε y^p against the matching description in Theorem 2.1.

    Returns (matches, which_statement).  This uses only the computed orbits —
    it does not consult the theorem's proof.
    """
    c, d, a, v = sol.c, sol.d, sol.alpha, sol.value
    vo = sol.value_orbit
    v_periodic = vo.preperiodic and vo.tail == 0
    per = vo.period if v_periodic else None
    if d == 2:
        ok = a in (v, -v) and v_periodic and per in (1, 2)
        return ok, "(1)"
    if d % 2 == 1:
        ok = a == v and v_periodic and per == 1
        return ok, "(2)"
    if c != -1:
        ok = a in (v, -v) and v_periodic and per == 1
        return ok, "(3)"
    # Statement (4) holds in the form f(α) = ε y^p, which is what the accepted version prints.
    # arXiv v1 prints α = ±ε y^p — the d = 2 relation, where N = 4 is even and f^4 fixes both
    # points of the 2-cycle {0, −1}.  For even d ≥ 4 the threshold N = 3 is odd, so f^3 *swaps*
    # the cycle: for f = x^4 − 1, f^3(0) = −1 and f^3(−1) = 0, and the v1 relation has no
    # solutions.  The conclusion (α preperiodic, ε y^p periodic) is the same either way.
    # ``statement_as_printed`` keeps the v1 reading so the difference stays visible in tests.
    ok = v in (0, -1) and v_periodic and per == 2 and a**d + c == v
    return ok, "(4)"


def statement_as_printed(sol: Solution) -> bool:
    """Theorem 2.1 (4) as printed in arXiv v1: α = ±ε y^p with ε y^p ∈ {0, −1} of exact period 2.
    Superseded in the accepted version; kept to document that the v1 reading fails for even
    d ≥ 4 (see ``statement``)."""
    vo = sol.value_orbit
    return sol.alpha in (sol.value, -sol.value) and sol.value in (0, -1) and vo.preperiodic and vo.tail == 0 and vo.period == 2


def search_theorem(d: int, c_values, alpha_bound: int, n: int | None = None) -> list[Solution]:
    """Every (c, α) with c in ``c_values`` (0 skipped), |α| ≤ alpha_bound and f^n(α) a prime-powered
    integer ε y^p (p | d).  ``n`` defaults to the theorem's threshold N."""
    n = threshold_iterate(d) if n is None else n
    found = []
    for c in c_values:
        if c == 0:
            continue
        for a in range(-alpha_bound, alpha_bound + 1):
            v = iterate(a, c, d, n)
            for eps, y, p in prime_power_forms(v, d):
                found.append(Solution(c, d, a, n, v, eps, y, p, orbit(a, c, d), orbit(v, c, d)))
    return found


def verify_theorem(d: int, c_values, alpha_bound: int) -> dict:
    """Run the search at the theorem's threshold and check every solution.  Returns a report."""
    sols = search_theorem(d, c_values, alpha_bound)
    bad = [s for s in sols if not (s.alpha_orbit.preperiodic and s.value_orbit.preperiodic
                                   and s.value_orbit.tail == 0 and statement(s)[0])]
    return {"d": d, "N": threshold_iterate(d), "solutions": sols, "violations": bad,
            "n_c": len([c for c in c_values if c != 0]), "alpha_bound": alpha_bound}


def search_near_misses(d: int, c_values, alpha_bound: int, n: int | None = None) -> list[Solution]:
    """Sharpness witnesses: f^n(α) = ε y^p with α NOT preperiodic, one iterate below the threshold.

    Remark 1.5 of the paper gives x² − 460 (f³(22) = 114²) and x^d − r^d (f²(r) = −r^d).
    """
    n = threshold_iterate(d) - 1 if n is None else n
    return [s for s in search_theorem(d, c_values, alpha_bound, n) if not s.alpha_orbit.preperiodic]


# --------------------------------------------------------------------------
# the small cases the paper settles by computer
# --------------------------------------------------------------------------

def mod8_obstruction(c: int, d: int = 2, n: int = 4) -> bool:
    """The Magma check in the proof of Theorem 2.1: for f = x² + c with c ∈ {1, 2},
    f⁴(α) ≡ ±y² (mod 8) has no solutions α, y ∈ Z/8Z — hence none in Z.

    Returns True when the congruence has no solutions (so Theorem 2.1 is vacuous for that c)."""
    residues = {(eps * y * y) % 8 for y in range(8) for eps in (1, -1)}
    for a in range(8):
        v = a
        for _ in range(n):
            v = (v**d + c) % 8
        if v in residues:
            return False
    return True
