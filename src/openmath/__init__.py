"""Numerical probes of results announced in the openai/math manuscript collection.

Each submodule corresponds to one result family in that collection and computes
the quantity appearing in the family's main theorem, so that the theorem's
predicted law can be compared against a finite numeration.

Nothing in this package attempts to verify a proof. It only evaluates the
objects the theorems are about.
"""

__all__ = [
    "sieves",
    "poisson_dirichlet",
    "joint_dickman",
    "jacobsthal",
    "prime_gaps",
    "egyptian",
    "pi_exponent",
    "quasi_riemann",
]
