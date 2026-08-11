#!/usr/bin/env python3
"""Exact symbolic derivations for quintic Hermite interpolation.

The script is deliberately both a derivation and a regression test.  All
polynomials and finite-difference coefficients are obtained from their defining
linear systems using exact SymPy arithmetic; the displayed formulas are then
generated from those objects.

Typical use::

    python derive_hermite.py --check
    python derive_hermite.py --emit-tex generated_formulas.tex
    python derive_hermite.py --check-tex generated_formulas.tex
"""

from __future__ import annotations

import argparse
import difflib
import itertools
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

import sympy as sp


Endpoint = int
DerivativeOrder = int
GridOffset = int


@dataclass(frozen=True)
class HermiteDerivation:
    """All exact one-dimensional objects used by the notes."""

    t: sp.Symbol
    psi: Tuple[sp.Expr, sp.Expr, sp.Expr]
    basis: Mapping[Tuple[Endpoint, DerivativeOrder], sp.Expr]
    jets: Mapping[Tuple[Endpoint, DerivativeOrder], sp.Symbol]
    interpolant: sp.Expr
    monomial_coefficients: Tuple[sp.Expr, ...]
    fd_offsets: Tuple[GridOffset, ...]
    fd_first: Tuple[sp.Rational, ...]
    fd_second: Tuple[sp.Rational, ...]
    grid_data: Mapping[GridOffset, sp.Symbol]
    value_only_interpolant: sp.Expr
    weights: Mapping[GridOffset, sp.Expr]
    weight_derivatives: Mapping[Tuple[GridOffset, DerivativeOrder], sp.Expr]
    moment_defects: Mapping[int, sp.Expr]


class Verifier:
    """Assertion helper that remains active even when Python uses ``-O``."""

    def __init__(self) -> None:
        self.count = 0

    def zero(self, expression: sp.Expr, description: str) -> None:
        self.count += 1
        if sp.expand(expression) != 0:
            simplified = sp.factor(expression)
            raise AssertionError(f"{description}: expected zero, got {simplified}")

    def equal(self, actual: object, expected: object, description: str) -> None:
        self.count += 1
        if isinstance(actual, sp.Basic) or isinstance(expected, sp.Basic):
            difference = sp.sympify(actual) - sp.sympify(expected)
            if sp.expand(difference) != 0:
                raise AssertionError(
                    f"{description}: expected {expected}, got {actual}"
                )
        elif actual != expected:
            raise AssertionError(f"{description}: expected {expected}, got {actual}")

    def true(self, condition: bool, description: str) -> None:
        self.count += 1
        if not condition:
            raise AssertionError(description)


def _solve_unique(equations: Sequence[sp.Expr], unknowns: Sequence[sp.Symbol]) -> Dict[sp.Symbol, sp.Expr]:
    """Solve a square exact system and require a unique complete solution."""

    solutions = sp.solve(equations, unknowns, dict=True, simplify=False)
    if len(solutions) != 1 or set(solutions[0]) != set(unknowns):
        raise RuntimeError("the defining linear system did not have one complete solution")
    return solutions[0]


def derive_left_cardinal_basis(t: sp.Symbol, derivative_order: int) -> sp.Expr:
    """Derive the left-end basis for one of value, slope, or curvature."""

    coefficients = sp.symbols(f"a0:{6}")
    polynomial = sum(coefficients[k] * t**k for k in range(6))
    equations = []
    for endpoint in (0, 1):
        for order in range(3):
            target = sp.Integer(endpoint == 0 and order == derivative_order)
            equations.append(
                sp.Eq(sp.diff(polynomial, t, order).subs(t, endpoint), target)
            )
    solution = _solve_unique(equations, coefficients)
    return sp.factor(polynomial.subs(solution))


def derive_fd_coefficients(
    offsets: Sequence[int], derivative_order: int
) -> Tuple[sp.Rational, ...]:
    """Derive a centered formula from polynomial moment equations.

    For five points, imposing moments zero through degree four gives the
    fourth-order first- and second-derivative formulas used in the paper.
    The returned coefficients approximate ``h**order * f**(order)``.
    """

    unknowns = sp.symbols(f"d0:{len(offsets)}")
    equations = []
    for degree in range(len(offsets)):
        target = sp.factorial(degree) if degree == derivative_order else 0
        equations.append(
            sp.Eq(
                sum(c * sp.Integer(j) ** degree for c, j in zip(unknowns, offsets)),
                target,
            )
        )
    solution = _solve_unique(equations, unknowns)
    return tuple(sp.Rational(solution[c]) for c in unknowns)


def derive_all() -> HermiteDerivation:
    """Construct every formula from its defining conditions."""

    t = sp.Symbol("t", real=True)
    psi = tuple(derive_left_cardinal_basis(t, order) for order in range(3))

    # Reflection changes the sign of odd derivatives.  The prefactor makes
    # every B_(e,r) cardinal with a *positive* r-th derivative at endpoint e.
    basis: Dict[Tuple[int, int], sp.Expr] = {}
    for endpoint in (0, 1):
        for order in range(3):
            argument = endpoint + (1 - 2 * endpoint) * t
            basis[endpoint, order] = sp.factor(
                (-1) ** (endpoint * order) * psi[order].subs(t, argument)
            )

    jets = {
        (endpoint, order): sp.Symbol(f"F_{{{endpoint},{order}}}")
        for endpoint in (0, 1)
        for order in range(3)
    }
    interpolant = sp.expand(
        sum(jets[key] * basis[key] for key in itertools.product((0, 1), range(3)))
    )
    monomial_coefficients = tuple(
        sp.expand(interpolant).coeff(t, power) for power in range(6)
    )

    fd_offsets = (-2, -1, 0, 1, 2)
    fd_first = derive_fd_coefficients(fd_offsets, 1)
    fd_second = derive_fd_coefficients(fd_offsets, 2)

    grid_data = {
        j: sp.Symbol(f"f_{{i{j:+d}}}" if j else "f_i") for j in range(-2, 4)
    }
    # Substitute h f' and h^2 f'' at each cell endpoint.  The powers of h
    # cancel because derive_fd_coefficients returns dimensionless stencils.
    value_only = grid_data[0] * basis[0, 0] + grid_data[1] * basis[1, 0]
    for coefficient, local_offset in zip(fd_first, fd_offsets):
        value_only += coefficient * grid_data[local_offset] * basis[0, 1]
        value_only += coefficient * grid_data[local_offset + 1] * basis[1, 1]
    for coefficient, local_offset in zip(fd_second, fd_offsets):
        value_only += coefficient * grid_data[local_offset] * basis[0, 2]
        value_only += coefficient * grid_data[local_offset + 1] * basis[1, 2]
    value_only = sp.expand(value_only)

    weights = {
        j: sp.factor(value_only.coeff(grid_data[j])) for j in range(-2, 4)
    }
    weight_derivatives = {
        (j, order): sp.factor(sp.diff(weights[j], t, order))
        for j in range(-2, 4)
        for order in range(3)
    }
    moment_defects = {
        degree: sp.factor(
            sum(weights[j] * sp.Integer(j) ** degree for j in range(-2, 4))
            - t**degree
        )
        for degree in range(11)
    }

    return HermiteDerivation(
        t=t,
        psi=psi,  # type: ignore[arg-type]
        basis=basis,
        jets=jets,
        interpolant=interpolant,
        monomial_coefficients=monomial_coefficients,
        fd_offsets=fd_offsets,
        fd_first=fd_first,
        fd_second=fd_second,
        grid_data=grid_data,
        value_only_interpolant=value_only,
        weights=weights,
        weight_derivatives=weight_derivatives,
        moment_defects=moment_defects,
    )


def fd_taylor_coefficient(
    coefficients: Sequence[sp.Expr], offsets: Sequence[int], derivative_order: int, degree: int
) -> sp.Expr:
    """Coefficient of ``h**(degree-order) f**(degree)`` in FD minus exact."""

    approximation = sum(
        coefficient * sp.Integer(offset) ** degree
        for coefficient, offset in zip(coefficients, offsets)
    ) / sp.factorial(degree)
    exact = sp.Integer(degree == derivative_order)
    return sp.factor(approximation - exact)


def _tensor_basis(
    derivation: HermiteDerivation,
    variables: Sequence[sp.Symbol],
    endpoints: Sequence[int],
    orders: Sequence[int],
) -> sp.Expr:
    return sp.prod(
        derivation.basis[endpoint, order].subs(derivation.t, variable)
        for variable, endpoint, order in zip(variables, endpoints, orders)
    )


def verify_one_dimensional(derivation: HermiteDerivation, check: Verifier) -> None:
    """Verify cardinality, exact Hermite reproduction, and the remainder."""

    t = derivation.t
    for endpoint, order in itertools.product((0, 1), range(3)):
        expected_reflection = (-1) ** (endpoint * order) * derivation.psi[
            order
        ].subs(t, endpoint + (1 - 2 * endpoint) * t)
        check.zero(
            derivation.basis[endpoint, order] - expected_reflection,
            f"oriented basis B({endpoint},{order})",
        )
        for evaluation_endpoint, evaluation_order in itertools.product(
            (0, 1), range(3)
        ):
            actual = sp.diff(
                derivation.basis[endpoint, order], t, evaluation_order
            ).subs(t, evaluation_endpoint)
            expected = sp.Integer(
                endpoint == evaluation_endpoint and order == evaluation_order
            )
            check.equal(
                actual,
                expected,
                "one-dimensional endpoint cardinality",
            )

    # The monomial coefficients were collected independently from H(t).
    reconstructed = sum(
        coefficient * t**power
        for power, coefficient in enumerate(derivation.monomial_coefficients)
    )
    check.zero(
        reconstructed - derivation.interpolant,
        "monomial and cardinal forms agree",
    )

    # Exact derivatives make the Hermite interpolant exact through degree five.
    for degree in range(6):
        polynomial = t**degree
        polynomial_interpolant = 0
        for endpoint, order in itertools.product((0, 1), range(3)):
            jet = sp.diff(polynomial, t, order).subs(t, endpoint)
            polynomial_interpolant += jet * derivation.basis[endpoint, order]
        check.zero(
            polynomial_interpolant - polynomial,
            f"exact Hermite reproduction of degree {degree}",
        )

    # For f(t)=t^6, f-H is exactly the nodal product in the remainder formula.
    degree_six = t**6
    degree_six_interpolant = sum(
        sp.diff(degree_six, t, order).subs(t, endpoint)
        * derivation.basis[endpoint, order]
        for endpoint, order in itertools.product((0, 1), range(3))
    )
    check.zero(
        degree_six - degree_six_interpolant - t**3 * (t - 1) ** 3,
        "degree-six Hermite remainder",
    )

    # Regression for the sign printed in SPHINCS_BSSN: using +psi_1(1-t)
    # for the right slope reverses the derivative there and fails for f(t)=t.
    wrong_linear = (
        derivation.psi[0].subs(t, 1 - t)
        + derivation.psi[1]
        + derivation.psi[1].subs(t, 1 - t)
    )
    check.true(
        sp.expand(wrong_linear - t) != 0,
        "the source's unoriented right-slope basis must fail for a line",
    )
    check.equal(
        sp.diff(wrong_linear, t).subs(t, 1),
        -1,
        "the source's right derivative has the reversed sign",
    )


def verify_finite_differences(derivation: HermiteDerivation, check: Verifier) -> None:
    """Verify defining moments and leading Taylor defects of both stencils."""

    for derivative_order, coefficients in (
        (1, derivation.fd_first),
        (2, derivation.fd_second),
    ):
        for degree in range(5):
            moment = sum(
                coefficient * sp.Integer(offset) ** degree
                for coefficient, offset in zip(coefficients, derivation.fd_offsets)
            )
            expected = sp.factorial(degree) if degree == derivative_order else 0
            check.equal(
                moment,
                expected,
                f"finite-difference moment q={derivative_order}, n={degree}",
            )

    check.equal(
        fd_taylor_coefficient(derivation.fd_first, derivation.fd_offsets, 1, 5),
        -sp.Rational(1, 30),
        "first-derivative leading Taylor defect",
    )
    check.equal(
        fd_taylor_coefficient(derivation.fd_first, derivation.fd_offsets, 1, 6),
        0,
        "first-derivative even Taylor defect",
    )
    check.equal(
        fd_taylor_coefficient(derivation.fd_second, derivation.fd_offsets, 2, 5),
        0,
        "second-derivative odd Taylor defect",
    )
    check.equal(
        fd_taylor_coefficient(derivation.fd_second, derivation.fd_offsets, 2, 6),
        -sp.Rational(1, 90),
        "second-derivative leading Taylor defect",
    )


def _endpoint_stencil(
    derivation: HermiteDerivation, endpoint: int, derivative_order: int
) -> Dict[int, sp.Expr]:
    t = derivation.t
    return {
        j: sp.diff(derivation.weights[j], t, derivative_order).subs(t, endpoint)
        for j in range(-2, 4)
    }


def verify_value_only_weights(derivation: HermiteDerivation, check: Verifier) -> None:
    """Verify the six weights, moment defects, midpoint, and C2 gluing."""

    t = derivation.t
    reconstructed = sum(
        derivation.weights[j] * derivation.grid_data[j] for j in range(-2, 4)
    )
    check.zero(
        reconstructed - derivation.value_only_interpolant,
        "collected six-point weights",
    )

    check.zero(
        sum(derivation.weights.values()) - 1,
        "partition of unity",
    )
    for order in range(3):
        for j in range(-2, 4):
            reflected = (-1) ** order * derivation.weight_derivatives[
                (1 - j, order)
            ].subs(t, 1 - t)
            check.zero(
                derivation.weight_derivatives[j, order] - reflected,
                f"weight reflection j={j}, q={order}",
            )

    # Reproduction through degree four, followed by the exact degree-five
    # defect.  Dividing the latter by 5! gives the leading local value error.
    for degree in range(5):
        check.zero(
            derivation.moment_defects[degree],
            f"six-point polynomial moment degree {degree}",
        )
    expected_degree_five_defect = (
        4 * t * (t - 1) * (2 * t - 1) * (3 * t**2 - 3 * t - 1)
    )
    check.zero(
        derivation.moment_defects[5] - expected_degree_five_defect,
        "degree-five moment defect",
    )
    leading_value_error = derivation.moment_defects[5] / sp.factorial(5)
    for derivative_order in range(3):
        check.true(
            sp.diff(leading_value_error, t, derivative_order) != 0,
            f"generic derivative q={derivative_order} has order h^(5-q)",
        )

    midpoint = sp.Rational(1, 2)
    midpoint_weights = (
        sp.Rational(3, 256),
        -sp.Rational(25, 256),
        sp.Rational(75, 128),
        sp.Rational(75, 128),
        -sp.Rational(25, 256),
        sp.Rational(3, 256),
    )
    for j, expected in zip(range(-2, 4), midpoint_weights):
        check.equal(
            derivation.weights[j].subs(t, midpoint),
            expected,
            f"midpoint weight j={j}",
        )
    check.equal(
        derivation.moment_defects[5].subs(t, midpoint),
        0,
        "midpoint fifth-degree superconvergence",
    )
    check.equal(
        derivation.moment_defects[6].subs(t, midpoint) / sp.factorial(6),
        sp.Rational(5, 1024),
        "midpoint leading sixth-derivative error",
    )
    check.equal(
        sp.diff(leading_value_error, t).subs(t, midpoint),
        sp.Rational(7, 240),
        "midpoint first-derivative leading error",
    )
    check.equal(
        sp.diff(leading_value_error, t, 2).subs(t, midpoint),
        0,
        "midpoint second derivative cancels the degree-five defect",
    )
    check.equal(
        sp.diff(
            derivation.moment_defects[6] / sp.factorial(6), t, 2
        ).subs(t, midpoint),
        -sp.Rational(259, 5760),
        "midpoint second-derivative leading sixth-degree error",
    )

    # At either endpoint, the first two t derivatives are exactly the centered
    # five-point stencils associated with that node.
    zero_value = {j: sp.Integer(j == 0) for j in range(-2, 4)}
    zero_first = {
        j: (derivation.fd_first[derivation.fd_offsets.index(j)] if j in derivation.fd_offsets else 0)
        for j in range(-2, 4)
    }
    zero_second = {
        j: (derivation.fd_second[derivation.fd_offsets.index(j)] if j in derivation.fd_offsets else 0)
        for j in range(-2, 4)
    }
    for order, expected_stencil in enumerate((zero_value, zero_first, zero_second)):
        actual_left = _endpoint_stencil(derivation, 0, order)
        for j in range(-2, 4):
            check.equal(
                actual_left[j],
                expected_stencil[j],
                f"left endpoint stencil q={order}, j={j}",
            )
            expected_right = expected_stencil.get(j - 1, 0)
            check.equal(
                _endpoint_stencil(derivation, 1, order)[j],
                expected_right,
                f"right endpoint shifted stencil q={order}, j={j}",
            )

    # Adjacent cells use the same nodal value and the same centered derivative
    # estimates.  Check equality of their value, slope, and curvature traces
    # directly in terms of independent grid samples f_{i-3},...,f_{i+3}.
    shared_data = {j: sp.Symbol(f"g_{{i{j:+d}}}") for j in range(-3, 4)}
    left_cell = sum(
        derivation.weights[j].subs(t, sp.Symbol("u")) * shared_data[j - 1]
        for j in range(-2, 4)
    )
    right_cell = sum(
        derivation.weights[j].subs(t, sp.Symbol("v")) * shared_data[j]
        for j in range(-2, 4)
    )
    u, v = sp.Symbol("u"), sp.Symbol("v")
    for order in range(3):
        left_trace = sp.diff(left_cell, u, order).subs(u, 1)
        right_trace = sp.diff(right_cell, v, order).subs(v, 0)
        check.zero(left_trace - right_trace, f"adjacent-cell C2 trace q={order}")


def verify_tensor_cardinality(
    derivation: HermiteDerivation, dimension: int, check: Verifier
) -> None:
    """Exhaustively verify all 36 or 216 corner-jet cardinal conditions."""

    endpoints = tuple(itertools.product((0, 1), repeat=dimension))
    orders = tuple(itertools.product(range(3), repeat=dimension))
    basis_indices = tuple(itertools.product(endpoints, orders))
    check.equal(
        len(basis_indices),
        6**dimension,
        f"{dimension}D tensor basis count",
    )
    check.equal(
        len(orders),
        3**dimension,
        f"{dimension}D derivative fields per corner",
    )

    # Precomputing the 36 univariate endpoint values makes the exhaustive
    # 216-by-216 3D check quick while still visiting every condition.
    endpoint_values = {
        (endpoint, order, evaluation_endpoint, evaluation_order): sp.diff(
            derivation.basis[endpoint, order], derivation.t, evaluation_order
        ).subs(derivation.t, evaluation_endpoint)
        for endpoint, order, evaluation_endpoint, evaluation_order in itertools.product(
            (0, 1), range(3), (0, 1), range(3)
        )
    }
    for (basis_endpoint, basis_order), (evaluation_endpoint, evaluation_order) in itertools.product(
        basis_indices, basis_indices
    ):
        actual = sp.prod(
            endpoint_values[e, r, ep, q]
            for e, r, ep, q in zip(
                basis_endpoint,
                basis_order,
                evaluation_endpoint,
                evaluation_order,
            )
        )
        expected = sp.Integer(
            basis_endpoint == evaluation_endpoint and basis_order == evaluation_order
        )
        check.equal(actual, expected, f"{dimension}D tensor cardinality")


def verify_tensor_structure(
    derivation: HermiteDerivation, dimension: int, check: Verifier
) -> None:
    """Verify degree, axis symmetry, mixed derivatives, and separability."""

    variables = sp.symbols(f"x0:{dimension}")
    endpoint_vectors = tuple(itertools.product((0, 1), repeat=dimension))
    order_vectors = tuple(itertools.product(range(3), repeat=dimension))

    univariate_degrees = {
        (endpoint, order): sp.Poly(
            derivation.basis[endpoint, order], derivation.t
        ).degree()
        for endpoint, order in itertools.product((0, 1), range(3))
    }

    for endpoints, orders in itertools.product(endpoint_vectors, order_vectors):
        factors = tuple(
            derivation.basis[endpoint, order].subs(derivation.t, variable)
            for variable, endpoint, order in zip(variables, endpoints, orders)
        )
        for endpoint, order in zip(endpoints, orders):
            check.equal(
                univariate_degrees[endpoint, order],
                5,
                f"{dimension}D tensor degree in each coordinate",
            )
        for axis in range(dimension - 1):
            # Check the factors after simultaneously exchanging the two axes.
            # This is equivalent to expanding both tensor products, but keeps
            # the 3D test linear in the number of basis functions.
            for destination_axis in (axis, axis + 1):
                source_axis = (
                    axis + 1 if destination_axis == axis else axis
                )
                exchanged_factor = factors[source_axis].subs(
                    variables[source_axis], variables[destination_axis]
                )
                expected_factor = derivation.basis[
                    endpoints[source_axis], orders[source_axis]
                ].subs(derivation.t, variables[destination_axis])
                check.zero(
                    exchanged_factor - expected_factor,
                    f"{dimension}D axis-exchange symmetry",
                )
        for first_axis in range(dimension):
            for second_axis in range(first_axis + 1, dimension):
                # Each coordinate occurs in one factor only, hence either
                # differentiation order produces these same exact factors.
                forward_factors = list(factors)
                reverse_factors = list(factors)
                forward_factors[first_axis] = sp.diff(
                    forward_factors[first_axis], variables[first_axis]
                )
                forward_factors[second_axis] = sp.diff(
                    forward_factors[second_axis], variables[second_axis]
                )
                reverse_factors[second_axis] = sp.diff(
                    reverse_factors[second_axis], variables[second_axis]
                )
                reverse_factors[first_axis] = sp.diff(
                    reverse_factors[first_axis], variables[first_axis]
                )
                check.true(
                    all(
                        forward_factor == reverse_factor
                        for forward_factor, reverse_factor in zip(
                            forward_factors, reverse_factors
                        )
                    ),
                    f"{dimension}D mixed-partial commutation",
                )

    # Every derivative of a separable value-only stencil weight is the product
    # of the corresponding 1D derivatives.  Visit all 36*9 and 216*27 cases.
    stencil_indices = tuple(itertools.product(range(-2, 4), repeat=dimension))
    derivative_vectors = tuple(itertools.product(range(3), repeat=dimension))
    factor_derivative_identity = {
        (axis, j, order): sp.expand(
            sp.diff(
                derivation.weights[j].subs(derivation.t, variables[axis]),
                variables[axis],
                order,
            )
            - derivation.weight_derivatives[j, order].subs(
                derivation.t, variables[axis]
            )
        )
        == 0
        for axis, j, order in itertools.product(
            range(dimension), range(-2, 4), range(3)
        )
    }
    for stencil_index, derivative_orders in itertools.product(
        stencil_indices, derivative_vectors
    ):
        # Differentiating a product whose factors use disjoint variables is a
        # factorwise operation.  Verify every factor for every tensor weight
        # and every derivative multi-index without expanding 3D polynomials.
        check.true(
            all(
                factor_derivative_identity[axis, j, order]
                for axis, (j, order) in enumerate(
                    zip(stencil_index, derivative_orders)
                )
            ),
            f"{dimension}D separated stencil derivative",
        )

    # Direct tensor summation and nested 1D interpolation are algebraically
    # identical.
    # In a nested interpolation, the coefficient of each independent datum is
    # obtained by following its one choice on every axis.  Exhaustively compare
    # those coefficients with the direct tensor-product coefficients.
    for index in stencil_indices:
        direct_coefficient_factors = tuple(
            derivation.weights[j].subs(derivation.t, variable)
            for j, variable in zip(index, variables)
        )
        nested_coefficient_factors = tuple(
            derivation.weights[index[axis]].subs(derivation.t, variables[axis])
            for axis in range(dimension)
        )
        check.true(
            direct_coefficient_factors == nested_coefficient_factors,
            f"{dimension}D direct/sequential coefficient identity",
        )

    # Traces on a shared face agree through normal derivative order two for
    # every normal and tangential cardinal basis factor.
    tangential_endpoints = tuple(itertools.product((0, 1), repeat=dimension - 1))
    tangential_orders = tuple(itertools.product(range(3), repeat=dimension - 1))
    for normal_basis_order, normal_trace_order, tangential_endpoint, tangential_order in itertools.product(
        range(3), range(3), tangential_endpoints, tangential_orders
    ):
        left_trace = (
            sp.diff(
                derivation.basis[1, normal_basis_order],
                derivation.t,
                normal_trace_order,
            ).subs(derivation.t, 1)
        )
        right_trace = (
            sp.diff(
                derivation.basis[0, normal_basis_order],
                derivation.t,
                normal_trace_order,
            ).subs(derivation.t, 0)
        )
        check.zero(
            left_trace - right_trace,
            f"{dimension}D shared-face C2 trace",
        )


def run_symbolic_checks(derivation: HermiteDerivation) -> int:
    """Run the complete exact regression suite and return assertion count."""

    check = Verifier()
    verify_one_dimensional(derivation, check)
    verify_finite_differences(derivation, check)
    verify_value_only_weights(derivation, check)
    for dimension in (2, 3):
        verify_tensor_cardinality(derivation, dimension, check)
        verify_tensor_structure(derivation, dimension, check)
    return check.count


def _latex(expression: sp.Expr) -> str:
    return sp.latex(sp.factor(expression))


def _align_block(rows: Iterable[Tuple[str, sp.Expr]]) -> list[str]:
    lines = [r"\begin{align*}"]
    rows = list(rows)
    for index, (left, right) in enumerate(rows):
        ending = r" \\" if index + 1 < len(rows) else ""
        lines.append(f"  {left} &= {_latex(right)}{ending}")
    lines.append(r"\end{align*}")
    return lines


def generate_tex(derivation: HermiteDerivation) -> str:
    """Return a deterministic, body-only LaTeX appendix fragment."""

    t = derivation.t
    lines = [
        "% Generated by derive_hermite.py; do not edit by hand.",
        "% Every displayed expression follows from an exact SymPy derivation.",
        "",
        r"\subsection{One-dimensional cardinal bases}",
        "",
        r"The left-end cardinal polynomials and the oriented endpoint bases are",
    ]
    lines.extend(
        _align_block(
            [(rf"\psi_{{{order}}}(t)", sp.expand(polynomial)) for order, polynomial in enumerate(derivation.psi)]
        )
    )
    lines.extend(
        [
            r"\[",
            r"  B_{e,r}(t)=(-1)^{er}\psi_r\!\left(e+(1-2e)t\right),",
            r"  \qquad e\in\{0,1\},\quad r\in\{0,1,2\}.",
            r"\]",
            r"Its normalized derivatives obey",
            r"\[",
            r"  B_{e,r}^{(q)}(t)=(-1)^{e(r+q)}",
            r"  \psi_r^{(q)}\!\left(e+(1-2e)t\right).",
            r"\]",
            "",
            r"\subsection{Derivatives of the cardinal polynomials}",
            "",
            r"\begingroup",
            r"\renewcommand{\arraystretch}{1.45}",
            r"\[\begin{array}{c|c|l}",
            r"r & q & \psi_r^{(q)}(t) \\ \hline",
        ]
    )
    for order, polynomial in enumerate(derivation.psi):
        for derivative_order in range(3):
            columns = [
                str(order),
                str(derivative_order),
                _latex(sp.expand(sp.diff(polynomial, t, derivative_order))),
            ]
            lines.append(" & ".join(columns) + r" \\")
        if order != 2:
            lines.append(r"\hline")
    lines.extend([r"\end{array}\]", r"\endgroup", ""])

    lines.extend(
        [
            r"\subsection{Monomial mapping}",
            "",
            r"With $F_{e,r}=h^r f^{(r)}(x_{i+e})$ and",
            r"$H(t)=\sum_{k=0}^{5}c_k t^k$, the coefficients are",
        ]
    )
    lines.extend(
        _align_block(
            [(rf"c_{{{power}}}", coefficient) for power, coefficient in enumerate(derivation.monomial_coefficients)]
        )
    )
    lines.append("")

    lines.extend(
        [
            r"\subsection{Uniform-grid six-point weights}",
            "",
            r"The derived fourth-order centered endpoint formulas are",
        ]
    )
    first_stencil = sum(
        coefficient * derivation.grid_data[offset]
        for coefficient, offset in zip(
            derivation.fd_first, derivation.fd_offsets
        )
    )
    second_stencil = sum(
        coefficient * derivation.grid_data[offset]
        for coefficient, offset in zip(
            derivation.fd_second, derivation.fd_offsets
        )
    )
    lines.extend(
        _align_block(
            [
                (r"h f'_i", first_stencil),
                (r"h^2 f''_i", second_stencil),
            ]
        )
    )
    lines.append(
        r"Collecting them in $H(t)=\sum_{j=-2}^{3}w_j(t)f_{i+j}$ gives"
    )
    lines.extend(
        _align_block(
            [(rf"w_{{{j}}}(t)", derivation.weights[j]) for j in range(-2, 4)]
        )
    )
    for derivative_order, ordinal, prime in (
        (1, "first", r"'"),
        (2, "second", r"''"),
    ):
        lines.extend(
            [
                "",
                rf"The explicit {ordinal} normalized derivatives are",
            ]
        )
        lines.extend(
            _align_block(
                [
                    (
                        rf"w{prime}_{{{j}}}(t)",
                        derivation.weight_derivatives[j, derivative_order],
                    )
                    for j in range(-2, 4)
                ]
            )
        )
    lines.extend(
        [
            "",
            r"Physical derivatives follow from",
            r"$\mathrm{d}^q H/\mathrm{d}x^q=h^{-q}\sum_j w_j^{(q)}(t)f_{i+j}$.",
            "",
        ]
    )
    return "\n".join(lines)


def write_tex(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def check_tex(path: Path, expected: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"generated TeX fragment does not exist: {path}")
    actual = path.read_text(encoding="utf-8")
    if actual != expected:
        difference = "".join(
            difflib.unified_diff(
                actual.splitlines(keepends=True),
                expected.splitlines(keepends=True),
                fromfile=str(path),
                tofile="fresh symbolic output",
            )
        )
        raise RuntimeError(f"generated TeX fragment is stale:\n{difference}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="derive and verify quintic Hermite interpolation formulas"
    )
    parser.add_argument(
        "--emit-tex",
        metavar="PATH",
        type=Path,
        help="write the deterministic generated LaTeX fragment",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="run all exact symbolic regression checks",
    )
    parser.add_argument(
        "--check-tex",
        metavar="PATH",
        type=Path,
        help="fail unless PATH exactly matches fresh generated LaTeX",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if not (arguments.emit_tex or arguments.check or arguments.check_tex):
        arguments.check = True

    derivation = derive_all()
    # Emission is only allowed after the same symbolic suite used by --check
    # passes, so the generated appendix can never bypass validation.
    assertion_count = run_symbolic_checks(derivation)

    generated = generate_tex(derivation)
    if arguments.emit_tex:
        write_tex(arguments.emit_tex, generated)
        print(f"wrote {arguments.emit_tex}")
    if arguments.check_tex:
        check_tex(arguments.check_tex, generated)
        print(f"verified {arguments.check_tex}")
    if arguments.check:
        print(f"all symbolic checks passed ({assertion_count} exact assertions)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
