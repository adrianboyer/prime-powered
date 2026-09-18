# prime-powered

**A computational companion I'm building to my own research paper.**

![python](https://img.shields.io/badge/python-3.10%2B-blue)
![license](https://img.shields.io/badge/license-MIT-lightgrey)
![status](https://img.shields.io/badge/status-in%20progress-yellow)

> A. Bhardwaj, A. Boyer-Paulet, W. Hindes, E. Qiu, A. Sun, *Prime-powered images and irreducible polynomials in dynamical semigroups*, [arXiv:2510.10310](https://arxiv.org/abs/2510.10310) — accepted at *Proc. Amer. Math. Soc. Ser. B*.

## Status

In progress. I'm building this myself, module by module, alongside coursework at Yale. Nothing below is finished yet — this file gets updated as pieces land.

## What it's for

The paper proves a classification theorem (Theorem 2.1) about what happens when you iterate `f(x) = x^d + c` and land on a prime power, and uses it to build irreducible polynomials in the semigroups `⟨x^d + c_1, …, x^d + c_s⟩`. I want a way to check the theorem computationally — brute-force search for every solution in a window and see whether it lines up with what the proof claims — rather than only trusting the proof by hand.

## Plan

- [ ] `dynamics.py` — iterate `f(x) = x^d + c` exactly; decide whether an integer point is preperiodic; search a window for solutions to `f^N(α) = ±y^p`
- [ ] `irreducibility.py` — check irreducibility of elements of the semigroup; verify the explicit irreducible families from Section 3
- [ ] tests for both, checked against the theorem's statements
- [ ] a short writeup of anything the search turns up

## Why

I want to understand every line of the proof well enough to defend it, not just cite it. Writing the brute-force check by hand is how I'm doing that.

MIT licensed.
