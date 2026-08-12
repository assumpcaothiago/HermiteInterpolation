"""Verify stable non-boosted LES Kerr formulas with SymPy.

The script follows Sec. 5.1 of T. Assumpcao, *A Hyperbolic Relaxation
Solver for Numerical Relativity* (2024). It assumes ``M > 0``,
``abs(a) < M``, and ``r > 0``. The LES radial coordinate double-covers
the Kerr exterior; ``r = 0`` is excluded because it is the puncture
representing the second asymptotic end.

Run ``python kerr_quasi_isotropic_check.py --show-components`` to display
the axis-regular Cartesian four-metric assembled by the checks.
"""

from __future__ import annotations

import argparse

import sympy as sp


def reduce_expression(expression: sp.Expr) -> sp.Expr:
    """Return a rationally simplified symbolic expression.

    The simplification proves rational identities only on their common
    nonsingular domain. Horizon and axis extensions are checked separately.

    :param expression: Expression to simplify.
    :return expression: Simplified expression.

    >>> u = sp.symbols("u")
    >>> reduce_expression((u**2 - 1) / (u - 1) - (u + 1))
    0
    """

    return sp.factor(sp.cancel(sp.together(expression)))


def require_zero(
    label: str, expression: sp.Expr | sp.MatrixBase
) -> None:
    """Raise an assertion error unless an expression or matrix is zero.

    :param label: Human-readable name of the check.
    :param expression: Scalar or matrix residual.
    :return: None.
    :raises AssertionError: If any simplified residual is nonzero.

    >>> require_zero("identity", sp.Integer(0))
    [PASS] identity
    """

    if isinstance(expression, sp.MatrixBase):
        residual = expression.applyfunc(reduce_expression)
        is_zero = residual == sp.zeros(*residual.shape)
    else:
        residual = reduce_expression(expression)
        is_zero = residual == 0
    if not is_zero:
        raise AssertionError(f"{label} failed:\n{residual}")
    print(f"[PASS] {label}")


def main(show_components: bool = False) -> None:
    """Construct the stable metric formulas and execute symbolic checks.

    :return: None.
    """

    # ------------------------------------------------------------------
    # 1. LES radial map, double covering, and stable horizon factors.
    # ------------------------------------------------------------------
    mass, radius = sp.symbols("M r", positive=True)
    spin = sp.symbols("a", real=True)
    # These are the recommended floating-point forms: the factorized radicand
    # avoids M**2 - a**2, while a**2/r_plus avoids the subtraction M - d.
    absolute_spin = sp.Abs(spin)
    horizon_gap = sp.sqrt(
        (mass - absolute_spin) * (mass + absolute_spin)
    )
    r_plus = mass + horizon_gap
    r_minus = spin**2 / r_plus
    r_minus_subtractive = mass - horizon_gap
    throat_radius = r_plus / 4

    r_boyer_map = (radius + throat_radius) ** 2 / radius
    p_factor = (radius - throat_radius) ** 2 / radius
    q_factor = p_factor + 2 * horizon_gap
    r_boyer_stable = r_plus + p_factor
    delta_stable = p_factor * q_factor
    radial_derivative = sp.diff(r_boyer_map, radius)

    require_zero("r_+ + r_- = 2M", r_plus + r_minus - 2 * mass)
    require_zero("r_+ r_- = a^2", r_plus * r_minus - spin**2)
    require_zero(
        "stable r_- agrees with M - d",
        r_minus - r_minus_subtractive,
    )
    require_zero("stable Boyer--Lindquist radius", r_boyer_stable - r_boyer_map)
    require_zero(
        "double covering R(r) = R(s^2/r)",
        r_boyer_map.subs(
            radius, throat_radius**2 / radius, simultaneous=True
        )
        - r_boyer_map,
    )
    require_zero("R - r_+ = P", r_boyer_map - r_plus - p_factor)
    require_zero("R - r_- = Q", r_boyer_map - r_minus - q_factor)
    require_zero(
        "dR/dr",
        radial_derivative
        - (radius - throat_radius)
        * (radius + throat_radius)
        / radius**2,
    )
    require_zero(
        "stable Delta = P Q",
        r_boyer_map**2
        - 2 * mass * r_boyer_map
        + spin**2
        - delta_stable,
    )
    require_zero(
        "radial metric transformation away from the throat",
        radial_derivative**2 / delta_stable
        - (radius + throat_radius) ** 2
        / (radius**3 * q_factor),
    )

    cosine_squared = sp.symbols("c2", nonnegative=True)
    sigma_radial = r_boyer_stable**2 + spin**2 * cosine_squared
    b_stable = (
        sigma_radial
        * (radius + throat_radius) ** 2
        / (radius**3 * q_factor)
    )
    sigma_at_throat = sigma_radial.subs(radius, throat_radius)
    require_zero(
        "finite subextremal horizon value of B",
        b_stable.subs(radius, throat_radius)
        - 4
        * sigma_at_throat
        / (throat_radius * (r_plus - r_minus)),
    )
    require_zero(
        "puncture scaling lim(r^4 B) = s^4",
        sp.limit(radius**4 * b_stable, radius, 0, dir="+")
        - throat_radius**4,
    )

    # ------------------------------------------------------------------
    # 2. Stable algebraic forms for A and B - Sigma/r^2.
    # ------------------------------------------------------------------
    boyer_radius, sine_squared = sp.symbols("R sin2", real=True)
    sigma_compact = boyer_radius**2 + spin**2 * (1 - sine_squared)
    delta_compact = boyer_radius**2 - 2 * mass * boyer_radius + spin**2
    a_boyer_lindquist = (
        (boyer_radius**2 + spin**2) ** 2
        - spin**2 * delta_compact * sine_squared
    )
    a_positive = (
        sigma_compact**2
        + spin**2
        * sine_squared
        * (sigma_compact + 2 * mass * boyer_radius)
    )
    require_zero(
        "cancellation-free form of A",
        a_boyer_lindquist - a_positive,
    )

    require_zero(
        "stable D = B - Sigma/r^2",
        b_stable
        - sigma_radial / radius**2
        - sigma_radial * r_minus / (radius**2 * q_factor),
    )

    # ------------------------------------------------------------------
    # 3. Inverse-coordinate Jacobian for (t,x,y,z) -> (t,r,theta,phi).
    # ------------------------------------------------------------------
    x_coord, y_coord, z_coord = sp.symbols("x y z", real=True)
    radius_squared = x_coord**2 + y_coord**2 + z_coord**2
    cylindrical_squared = x_coord**2 + y_coord**2
    radius_xyz = sp.sqrt(radius_squared)
    cylindrical_radius = sp.sqrt(cylindrical_squared)

    jacobian = sp.Matrix(
        [
            [1, 0, 0, 0],
            [
                0,
                x_coord / radius_xyz,
                y_coord / radius_xyz,
                z_coord / radius_xyz,
            ],
            [
                0,
                x_coord * z_coord / (radius_squared * cylindrical_radius),
                y_coord * z_coord / (radius_squared * cylindrical_radius),
                -cylindrical_radius / radius_squared,
            ],
            [
                0,
                -y_coord / cylindrical_squared,
                x_coord / cylindrical_squared,
                0,
            ],
        ]
    )

    theta_xyz = sp.acos(z_coord / radius_xyz)
    phi_xyz = sp.atan2(y_coord, x_coord)
    spatial_coordinates = (x_coord, y_coord, z_coord)
    differentiated_jacobian = sp.eye(4)
    differentiated_jacobian[1, 1:4] = sp.Matrix(
        1, 3, [sp.diff(radius_xyz, value) for value in spatial_coordinates]
    )
    differentiated_jacobian[2, 1:4] = sp.Matrix(
        1, 3, [sp.diff(theta_xyz, value) for value in spatial_coordinates]
    )
    differentiated_jacobian[3, 1:4] = sp.Matrix(
        1, 3, [sp.diff(phi_xyz, value) for value in spatial_coordinates]
    )
    require_zero(
        "inverse-coordinate Jacobian off the axis",
        differentiated_jacobian - jacobian,
    )

    # ------------------------------------------------------------------
    # 4. Exact basis change and removal of all polar-axis denominators.
    # ------------------------------------------------------------------
    g_time_time, h_function = sp.symbols("G_tt H")
    b_function, l_function, c_function = sp.symbols("B L C")
    d_function = b_function - l_function

    # The identity A - Sigma^2 implies lambda_phi = L + C varpi^2.
    spherical_metric = sp.Matrix(
        [
            [g_time_time, 0, 0, -h_function * cylindrical_squared],
            [0, b_function, 0, 0],
            [0, 0, l_function * radius_squared, 0],
            [
                -h_function * cylindrical_squared,
                0,
                0,
                (l_function + c_function * cylindrical_squared)
                * cylindrical_squared,
            ],
        ]
    )
    metric_from_jacobian = jacobian.T * spherical_metric * jacobian

    n_vector = sp.Matrix([x_coord, y_coord, z_coord]) / radius_xyz
    m_vector = sp.Matrix([-y_coord, x_coord, 0])
    gamma_regular = (
        l_function * sp.eye(3)
        + d_function * n_vector * n_vector.T
        + c_function * m_vector * m_vector.T
    )
    metric_regular = sp.zeros(4)
    metric_regular[0, 0] = g_time_time
    metric_regular[0, 1] = metric_regular[1, 0] = h_function * y_coord
    metric_regular[0, 2] = metric_regular[2, 0] = -h_function * x_coord
    metric_regular[1:4, 1:4] = gamma_regular
    require_zero(
        "axis-regular Cartesian four-metric from J^T g J",
        metric_from_jacobian - metric_regular,
    )
    generic_four_metric_determinant = (
        b_function
        * l_function
        * (
            g_time_time
            * (l_function + c_function * cylindrical_squared)
            - h_function**2 * cylindrical_squared
        )
    )
    require_zero(
        "axis-regular four-metric determinant",
        metric_regular.det() - generic_four_metric_determinant,
    )
    require_zero(
        "axis-regular spatial determinant",
        gamma_regular.det()
        - b_function
        * l_function
        * (l_function + c_function * cylindrical_squared),
    )

    positive_axis_radius = sp.symbols("z_axis", positive=True)
    gamma_on_axis = gamma_regular.subs(
        {
            x_coord: 0,
            y_coord: 0,
            z_coord: positive_axis_radius,
        },
        simultaneous=True,
    )
    require_zero(
        "direct positive polar-axis values",
        gamma_on_axis - sp.diag(l_function, l_function, b_function),
    )
    gamma_on_negative_axis = gamma_regular.subs(
        {
            x_coord: 0,
            y_coord: 0,
            z_coord: -positive_axis_radius,
        },
        simultaneous=True,
    )
    require_zero(
        "direct negative polar-axis values",
        gamma_on_negative_axis
        - sp.diag(l_function, l_function, b_function),
    )

    # ------------------------------------------------------------------
    # 5. Fully substituted stable Cartesian metric and ADM identities.
    # ------------------------------------------------------------------
    p_xyz = (radius_xyz - throat_radius) ** 2 / radius_xyz
    q_xyz = p_xyz + 2 * horizon_gap
    r_boyer_xyz = r_plus + p_xyz
    delta_xyz = p_xyz * q_xyz
    sigma_xyz = r_boyer_xyz**2 + spin**2 * z_coord**2 / radius_squared
    a_xyz = (
        sigma_xyz**2
        + spin**2
        * cylindrical_squared
        / radius_squared
        * (sigma_xyz + 2 * mass * r_boyer_xyz)
    )
    b_xyz = (
        sigma_xyz
        * (radius_xyz + throat_radius) ** 2
        / (radius_xyz**3 * q_xyz)
    )
    l_xyz = sigma_xyz / radius_squared
    d_xyz = sigma_xyz * r_minus / (radius_squared * q_xyz)
    c_xyz = (
        spin**2
        * (sigma_xyz + 2 * mass * r_boyer_xyz)
        / (sigma_xyz * radius_xyz**4)
    )
    h_xyz = 2 * mass * spin * r_boyer_xyz / (sigma_xyz * radius_squared)

    gamma_xyz = (
        l_xyz * sp.eye(3)
        + d_xyz * n_vector * n_vector.T
        + c_xyz * m_vector * m_vector.T
    )
    beta_up_xyz = sp.Matrix(
        [
            2 * mass * spin * r_boyer_xyz * y_coord / a_xyz,
            -2 * mass * spin * r_boyer_xyz * x_coord / a_xyz,
            0,
        ]
    )
    beta_down_xyz = sp.Matrix(
        [h_xyz * y_coord, -h_xyz * x_coord, 0]
    )
    alpha_squared_xyz = delta_xyz * sigma_xyz / a_xyz
    g_time_time_xyz = -(1 - 2 * mass * r_boyer_xyz / sigma_xyz)

    require_zero(
        "physical azimuthal spatial eigenvalue",
        l_xyz
        + c_xyz * cylindrical_squared
        - a_xyz / (sigma_xyz * radius_squared),
    )
    require_zero(
        "physical Cartesian spatial determinant",
        b_xyz * a_xyz / radius_xyz**4
        - (
            b_xyz
            * l_xyz
            * (l_xyz + c_xyz * cylindrical_squared)
        ),
    )
    # The generic determinant checked above is B*L times the bracket below.
    # Since L=Sigma/r^2, this residual proves
    # det(g)=-Delta*Sigma*B/r^4 without expanding a large 4-by-4 determinant
    # after all physical substitutions.
    require_zero(
        "fully substituted four-metric determinant factor",
        g_time_time_xyz
        * (l_xyz + c_xyz * cylindrical_squared)
        - h_xyz**2 * cylindrical_squared
        + delta_xyz / radius_squared,
    )
    require_zero(
        "fully substituted gamma_ij beta^j = beta_i",
        gamma_xyz * beta_up_xyz - beta_down_xyz,
    )
    require_zero(
        "fully substituted -alpha^2 + beta_i beta^i = g_tt",
        -alpha_squared_xyz
        + (beta_down_xyz.T * beta_up_xyz)[0]
        - g_time_time_xyz,
    )

    four_metric_determinant_radial = (
        -delta_stable * sigma_radial * b_stable / radius**4
    )
    require_zero(
        "four-metric determinant vanishes at the throat",
        four_metric_determinant_radial.subs(radius, throat_radius),
    )

    explicit_metric = sp.zeros(4)
    explicit_metric[0, 0] = g_time_time_xyz
    explicit_metric[0, 1] = explicit_metric[1, 0] = h_xyz * y_coord
    explicit_metric[0, 2] = explicit_metric[2, 0] = -h_xyz * x_coord
    explicit_metric[1:4, 1:4] = gamma_xyz

    if show_components:
        print("\nAxis-regular Cartesian four-metric g_(mu nu):")
        sp.pprint(explicit_metric)


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    :return: Parsed command-line arguments.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--show-components",
        action="store_true",
        help="print the assembled axis-regular Cartesian four-metric",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_arguments()
    main(show_components=arguments.show_components)
