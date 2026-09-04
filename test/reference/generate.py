#!/usr/bin/env python3
"""Generate deterministic multiprecision Cauchy-transform references."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import mpmath as mp


MPMATH_VERSION = "1.3.0"
WORKING_DIGITS = 400
OUTPUT_DIGITS = 60


@dataclass(frozen=True)
class Model:
    name: str
    degree: int
    knots: tuple
    coefficients: tuple
    samples: tuple = ()


def decimal(text: str):
    return mp.mpf(text)


def polynomial_value(coefficients, argument):
    result = mp.mpf("0")
    for coefficient in reversed(coefficients):
        result = result * argument + coefficient
    return result


def polynomial_derivative(coefficients, argument):
    result = mp.mpf("0")
    for power in range(len(coefficients) - 1, 0, -1):
        result = result * argument + power * coefficients[power]
    return result


def natural_cubic_coefficients(knots, values):
    """Return natural-spline coefficients in u=(E-left)/width."""
    count = len(knots)
    diagonal = []
    lower = []
    upper = []
    right_hand_side = []
    for point in range(1, count - 1):
        left_width = knots[point] - knots[point - 1]
        right_width = knots[point + 1] - knots[point]
        lower.append(left_width)
        diagonal.append(2 * (left_width + right_width))
        upper.append(right_width)
        right_hand_side.append(
            6
            * (
                (values[point + 1] - values[point]) / right_width
                - (values[point] - values[point - 1]) / left_width
            )
        )

    # Deterministic Thomas solve for the interior second derivatives.
    for row in range(1, len(diagonal)):
        factor = lower[row] / diagonal[row - 1]
        diagonal[row] -= factor * upper[row - 1]
        right_hand_side[row] -= factor * right_hand_side[row - 1]
    interior = [mp.mpf("0")] * len(diagonal)
    if interior:
        interior[-1] = right_hand_side[-1] / diagonal[-1]
        for row in range(len(interior) - 2, -1, -1):
            interior[row] = (
                right_hand_side[row] - upper[row] * interior[row + 1]
            ) / diagonal[row]
    second_derivatives = [mp.mpf("0"), *interior, mp.mpf("0")]

    intervals = []
    for interval in range(count - 1):
        width = knots[interval + 1] - knots[interval]
        delta = values[interval + 1] - values[interval]
        second_left = second_derivatives[interval]
        second_right = second_derivatives[interval + 1]
        intervals.append(
            (
                values[interval],
                delta - width**2 * (2 * second_left + second_right) / 6,
                width**2 * second_left / 2,
                width**2 * (second_right - second_left) / 6,
            )
        )
    return tuple(intervals)


def akima_coefficients(knots, values):
    """Return Akima-spline coefficients in u=(E-left)/width."""
    slopes = [
        (values[index + 1] - values[index]) / (knots[index + 1] - knots[index])
        for index in range(len(knots) - 1)
    ]
    extended = [
        3 * slopes[0] - 2 * slopes[1],
        2 * slopes[0] - slopes[1],
        *slopes,
        2 * slopes[-1] - slopes[-2],
        3 * slopes[-1] - 2 * slopes[-2],
    ]
    derivatives = []
    for point in range(len(knots)):
        left_weight = abs(extended[point + 3] - extended[point + 2])
        right_weight = abs(extended[point + 1] - extended[point])
        if left_weight + right_weight == 0:
            derivatives.append((extended[point + 1] + extended[point + 2]) / 2)
        else:
            derivatives.append(
                (
                    left_weight * extended[point + 1]
                    + right_weight * extended[point + 2]
                )
                / (left_weight + right_weight)
            )

    intervals = []
    for interval in range(len(knots) - 1):
        width = knots[interval + 1] - knots[interval]
        delta = values[interval + 1] - values[interval]
        intervals.append(
            (
                values[interval],
                width * derivatives[interval],
                3 * delta
                - width * (2 * derivatives[interval] + derivatives[interval + 1]),
                -2 * delta
                + width * (derivatives[interval] + derivatives[interval + 1]),
            )
        )
    return tuple(intervals)


def make_models():
    quadratic_knots = tuple(
        map(decimal, ("-1.875", "-1.125", "0.25", "1.875"))
    )
    quadratic_coefficients = tuple(
        tuple(map(decimal, row))
        for row in (
            ("0.25", "0.625", "-0.125"),
            ("0.75", "-0.375", "0.125"),
            ("0.5", "0.5625", "-0.1875"),
        )
    )

    cubic_samples = tuple(
        (decimal(energy), decimal(value))
        for energy, value in (
            ("-2.0", "0.125"),
            ("-1.375", "0.5625"),
            ("-0.25", "0.9375"),
            ("0.625", "0.6875"),
            ("1.25", "0.3125"),
            ("2.125", "0.0625"),
        )
    )
    cubic_knots = tuple(point[0] for point in cubic_samples)
    cubic_values = tuple(point[1] for point in cubic_samples)
    cubic_coefficients = natural_cubic_coefficients(cubic_knots, cubic_values)

    legendre_knots = tuple(map(decimal, ("-1.0", "1.0")))
    # P3(E) = (5 E^3 - 3 E) / 2 in u=(E+1)/2 coordinates.
    legendre_coefficients = (tuple(map(decimal, ("-1.0", "12.0", "-30.0", "20.0"))),)

    cancellation_knots = tuple(map(decimal, ("0.0", "1.0")))
    cancellation_coefficients = (tuple(map(decimal, ("1e100", "-2e100", "1.0"))),)

    return (
        Model("quadratic", 2, quadratic_knots, quadratic_coefficients),
        Model(
            "cubic_cspline",
            3,
            cubic_knots,
            cubic_coefficients,
            cubic_samples,
        ),
        Model("legendre3", 3, legendre_knots, legendre_coefficients),
        Model("ill_conditioned", 2, cancellation_knots, cancellation_coefficients),
    )


def weighted_coefficients(model, interval, power):
    left = model.knots[interval]
    width = model.knots[interval + 1] - left
    monomial = tuple(
        mp.binomial(power, local_power)
        * left ** (power - local_power)
        * width**local_power
        for local_power in range(power + 1)
    )
    source = model.coefficients[interval]
    product = [mp.mpf("0")] * (len(source) + power)
    for first, coefficient in enumerate(source):
        for second, factor in enumerate(monomial):
            product[first + second] += coefficient * factor
    return tuple(product)


def interval_transform(coefficients, root, complex_argument):
    if complex_argument:
        base = mp.log(root) - mp.log(root - 1)
    else:
        base = mp.log(abs(root)) - mp.log(abs(root - 1))
    result = coefficients[0] * base
    monomial_transform = base
    for power in range(1, len(coefficients)):
        monomial_transform = root * monomial_transform - mp.mpf(1) / power
        result += coefficients[power] * monomial_transform
    return result


def cauchy_transform(model, power, argument):
    result = mp.mpc(0)
    for interval in range(len(model.knots) - 1):
        left = model.knots[interval]
        width = model.knots[interval + 1] - left
        root = (argument - left) / width
        result += interval_transform(
            weighted_coefficients(model, interval, power), root, True
        )
    return result


def model_value(model, argument):
    if argument == model.knots[-1]:
        interval = len(model.knots) - 2
    else:
        interval = next(
            index
            for index in range(len(model.knots) - 1)
            if model.knots[index] <= argument < model.knots[index + 1]
        )
    left = model.knots[interval]
    width = model.knots[interval + 1] - left
    return polynomial_value(model.coefficients[interval], (argument - left) / width)


def divided_remainder_integral(coefficients, root, value_at_root):
    adjusted = list(coefficients)
    adjusted[0] -= value_at_root
    if len(adjusted) == 1:
        remainder = adjusted[0]
        quotient = []
    else:
        quotient = [mp.mpf("0")] * (len(adjusted) - 1)
        quotient[-1] = adjusted[-1]
        for power in range(len(quotient) - 1, 0, -1):
            quotient[power - 1] = adjusted[power] + root * quotient[power]
        remainder = adjusted[0] + root * quotient[0]
    scale = max([mp.mpf("1"), *map(abs, adjusted)])
    if abs(remainder) > decimal("1e-100") * scale:
        raise ArithmeticError("PV subtraction did not cancel at the singularity")
    return -sum(
        coefficient / (power + 1) for power, coefficient in enumerate(quotient)
    )


def principal_value_transform(model, power, argument):
    inside = model.knots[0] < argument < model.knots[-1]
    if not inside:
        result = mp.mpf("0")
        for interval in range(len(model.knots) - 1):
            left = model.knots[interval]
            width = model.knots[interval + 1] - left
            root = (argument - left) / width
            result += interval_transform(
                weighted_coefficients(model, interval, power), root, False
            )
        return result

    value_at_argument = argument**power * model_value(model, argument)
    result = value_at_argument * (
        mp.log(abs(argument - model.knots[0]))
        - mp.log(abs(argument - model.knots[-1]))
    )
    for interval in range(len(model.knots) - 1):
        left = model.knots[interval]
        right = model.knots[interval + 1]
        width = right - left
        root = (argument - left) / width
        coefficients = weighted_coefficients(model, interval, power)
        if left <= argument <= right:
            result += divided_remainder_integral(
                coefficients, root, value_at_argument
            )
        else:
            adjusted = list(coefficients)
            adjusted[0] -= value_at_argument
            result += interval_transform(adjusted, root, False)
    return result


def moment(model, power):
    terms = []
    for interval in range(len(model.knots) - 1):
        width = model.knots[interval + 1] - model.knots[interval]
        coefficients = weighted_coefficients(model, interval, power)
        terms.extend(
            width * coefficient / (local_power + 1)
            for local_power, coefficient in enumerate(coefficients)
        )
    result = mp.fsum(terms)
    if abs(result) <= 8 * mp.eps * mp.fsum(map(abs, terms)):
        return mp.mpf("0")
    return result


def numerical_complex_transform(model, power, argument):
    result = mp.mpc(0)
    for interval in range(len(model.knots) - 1):
        left = model.knots[interval]
        width = model.knots[interval + 1] - left
        coefficients = weighted_coefficients(model, interval, power)
        root = (argument - left) / width
        result += mp.quad(
            lambda local: polynomial_value(coefficients, local) / (root - local),
            [0, 1],
        )
    return result


def numerical_principal_value(model, power, argument):
    inside = model.knots[0] < argument < model.knots[-1]
    value_at_argument = (
        argument**power * model_value(model, argument) if inside else mp.mpf("0")
    )
    result = mp.mpf("0")
    for interval in range(len(model.knots) - 1):
        left = model.knots[interval]
        right = model.knots[interval + 1]
        width = right - left
        root = (argument - left) / width
        coefficients = weighted_coefficients(model, interval, power)

        def integrand(local):
            if local == root:
                return -polynomial_derivative(coefficients, root)
            return (
                polynomial_value(coefficients, local) - value_at_argument
            ) / (root - local)

        points = [mp.mpf("0"), mp.mpf("1")]
        if 0 < root < 1:
            points.insert(1, root)
        result += mp.quad(integrand, points)
    if inside:
        result += value_at_argument * (
            mp.log(abs(argument - model.knots[0]))
            - mp.log(abs(argument - model.knots[-1]))
        )
    return result


def assert_close(actual, expected, description):
    if not mp.isfinite(actual) or not mp.isfinite(expected):
        raise ArithmeticError(
            f"multiprecision cross-check produced a non-finite value for {description}"
        )
    tolerance = decimal("1e-85") * max(decimal("1e-300"), abs(actual), abs(expected))
    if abs(actual - expected) > tolerance:
        raise ArithmeticError(
            f"multiprecision cross-check failed for {description}: "
            f"error={mp.nstr(abs(actual - expected), 12)}"
        )


def validate_models(models):
    for model in models:
        for interval, coefficients in enumerate(model.coefficients):
            if len(coefficients) != model.degree + 1:
                raise ArithmeticError(f"wrong degree for {model.name}")
            if interval + 1 < len(model.coefficients):
                assert_close(
                    polynomial_value(coefficients, 1),
                    model.coefficients[interval + 1][0],
                    f"{model.name} continuity at knot {interval + 1}",
                )
        for energy, value in model.samples:
            assert_close(
                model_value(model, energy),
                value,
                f"{model.name} interpolation at {energy}",
            )


def format_number(value):
    if value == 0:
        return "0.0"
    return mp.nstr(value, OUTPUT_DIGITS, strip_zeros=False)


def make_table(columns, rows):
    lines = [
        "# generated by test/reference/generate.py; do not edit",
        f"# working_decimal_digits: {WORKING_DIGITS}",
        f"# output_significant_digits: {OUTPUT_DIGITS}",
        "# columns: " + "\t".join(columns),
    ]
    lines.extend("\t".join(row) for row in rows)
    return ("\n".join(lines) + "\n").encode("ascii")


def make_plain_rows(rows):
    return ("\n".join(" ".join(row) for row in rows) + "\n").encode("ascii")


def build_outputs():
    models = make_models()
    by_name = {model.name: model for model in models}
    validate_models(models)

    kk_samples = tuple(
        (decimal(energy), decimal(value))
        for energy, value in (
            ("-5", "0.5"),
            ("-3", "-2"),
            ("-1.25", "1.5"),
            ("-0.4", "0.25"),
            ("0.4", "3"),
            ("1.25", "-1"),
            ("3", "2.5"),
            ("5", "-0.75"),
        )
    )
    kk_knots = tuple(point[0] for point in kk_samples)
    kk_values = tuple(point[1] for point in kk_samples)
    kk_models = (
        Model("cspline", 3, kk_knots, natural_cubic_coefficients(kk_knots, kk_values), kk_samples),
        Model("akima", 3, kk_knots, akima_coefficients(kk_knots, kk_values), kk_samples),
    )
    validate_models(kk_models)

    complex_cases = (
        ("quadratic_upper_knot", "quadratic", 0, "-1.125", "0.125"),
        ("quadratic_lower_knot", "quadratic", 0, "-1.125", "-0.125"),
        ("quadratic_power2", "quadratic", 2, "0.3125", "0.0625"),
        ("quadratic_power3_lower_outside", "quadratic", 3, "2.5", "-0.375"),
        ("cubic_upper_knot", "cubic_cspline", 0, "-0.25", "0.09375"),
        ("cubic_lower_knot", "cubic_cspline", 0, "-0.25", "-0.09375"),
        ("cubic_power1", "cubic_cspline", 1, "0.75", "0.25"),
        ("cubic_power2_lower", "cubic_cspline", 2, "-1.5", "-0.1875"),
        ("legendre3_far", "legendre3", 0, "1e16", "2e15"),
        ("legendre3_threshold_near", "legendre3", 0, "2.998", "0.125"),
        ("legendre3_threshold_far", "legendre3", 0, "2.9981", "0.125"),
        ("ill_conditioned_far", "ill_conditioned", 0, "1e101", "2e100"),
    )
    real_cases = (
        ("quadratic_pv_knot_left", "quadratic", 0, "-1.125"),
        ("quadratic_pv_knot_right_power2", "quadratic", 2, "0.25"),
        ("quadratic_outside_below_power1", "quadratic", 1, "-2.75"),
        ("quadratic_outside_above_power3", "quadratic", 3, "2.5"),
        ("cubic_pv_knot_left", "cubic_cspline", 0, "-1.375"),
        ("cubic_pv_knot_right_power3", "cubic_cspline", 3, "0.625"),
        ("cubic_pv_interior_power1", "cubic_cspline", 1, "0.875"),
        ("cubic_outside_below_power2", "cubic_cspline", 2, "-2.75"),
        ("cubic_outside_above", "cubic_cspline", 0, "3.0"),
        ("legendre3_far_real", "legendre3", 0, "1e16"),
        ("ill_conditioned_far_real", "ill_conditioned", 0, "1e101"),
    )

    polynomial_rows = []
    for model in models:
        for interval, coefficients in enumerate(model.coefficients):
            padded = (*coefficients, *([mp.mpf("0")] * (4 - len(coefficients))))
            polynomial_rows.append(
                (
                    model.name,
                    str(interval),
                    format_number(model.knots[interval]),
                    format_number(model.knots[interval + 1]),
                    str(model.degree),
                    *(format_number(value) for value in padded),
                )
            )

    moment_rows = []
    for model in models:
        for power in range(4):
            moment_rows.append(
                (model.name, str(power), format_number(moment(model, power)))
            )

    transform_rows = []
    for case, model_name, power, real_part, imaginary_part in complex_cases:
        model = by_name[model_name]
        argument = mp.mpc(decimal(real_part), decimal(imaginary_part))
        result = cauchy_transform(model, power, argument)
        assert_close(
            result,
            numerical_complex_transform(model, power, argument),
            case,
        )
        transform_rows.append(
            (
                case,
                model_name,
                str(power),
                "complex",
                format_number(argument.real),
                format_number(argument.imag),
                format_number(result.real),
                format_number(result.imag),
            )
        )
    for case, model_name, power, real_part in real_cases:
        model = by_name[model_name]
        argument = decimal(real_part)
        result = principal_value_transform(model, power, argument)
        assert_close(
            result,
            numerical_principal_value(model, power, argument),
            case,
        )
        transform_rows.append(
            (
                case,
                model_name,
                str(power),
                "real_pv",
                format_number(argument),
                "0.0",
                format_number(result),
                "0.0",
            )
        )

    cubic = by_name["cubic_cspline"]
    omega = decimal("0.375")
    sigma_real = decimal("-0.125")
    sigma_imaginary = decimal("-0.25")
    dmft_argument = mp.mpc(omega - sigma_real, -sigma_imaginary)
    dmft_transform = cauchy_transform(cubic, 0, dmft_argument)
    assert_close(
        dmft_transform,
        numerical_complex_transform(cubic, 0, dmft_argument),
        "dmft_basic",
    )
    mapped = -dmft_transform / mp.pi
    dmft_rows = [
        (
            "dmft_basic",
            cubic.name,
            format_number(omega),
            format_number(sigma_real),
            format_number(sigma_imaginary),
            format_number(dmft_argument.real),
            format_number(dmft_argument.imag),
            format_number(dmft_transform.real),
            format_number(dmft_transform.imag),
            format_number(mapped.real),
            format_number(mapped.imag),
        )
    ]

    kk_rows = []
    for model in kk_models:
        for probe in map(decimal, ("-4", "-2", "0", "2", "4")):
            canonical = principal_value_transform(model, 0, probe)
            assert_close(canonical, numerical_principal_value(model, 0, probe), f"kk_{model.name}_{probe}")
            kk_rows.append((model.name, format_number(probe), format_number(-canonical / mp.pi)))

    return {
        "polynomials.tsv": make_table(
            ("model", "interval", "left", "right", "degree", "c0", "c1", "c2", "c3"),
            polynomial_rows,
        ),
        "moments.tsv": make_table(("model", "power", "moment"), moment_rows),
        "cauchy_transforms.tsv": make_table(
            (
                "case",
                "model",
                "power",
                "argument_kind",
                "z_re",
                "z_im",
                "expected_re",
                "expected_im",
            ),
            transform_rows,
        ),
        "dmft_mapping.tsv": make_table(
            (
                "case",
                "model",
                "omega",
                "sigma_re",
                "sigma_im",
                "z_re",
                "z_im",
                "h_re",
                "h_im",
                "mapped_re",
                "mapped_im",
            ),
            dmft_rows,
        ),
        "cubic_cspline_dos.dat": make_plain_rows(
            (format_number(energy), format_number(value))
            for energy, value in cubic.samples
        ),
        "dmft_resigma.dat": make_plain_rows(
            [(format_number(omega), format_number(sigma_real))]
        ),
        "dmft_imsigma.dat": make_plain_rows(
            [(format_number(omega), format_number(sigma_imaginary))]
        ),
        "dmft_reaw.ref": make_plain_rows(
            [(format_number(omega), format_number(mapped.real))]
        ),
        "dmft_imaw.ref": make_plain_rows(
            [(format_number(omega), format_number(mapped.imag))]
        ),
        "kk_nonlinear_input.dat": make_plain_rows(
            (format_number(energy), format_number(value)) for energy, value in kk_samples
        ),
        "kk_nonlinear.tsv": make_table(("method", "energy", "expected_minus_pv_over_pi"), kk_rows),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare generated bytes with the committed files without writing",
    )
    arguments = parser.parse_args()

    if mp.__version__ != MPMATH_VERSION:
        raise SystemExit(
            f"mpmath {MPMATH_VERSION} is required; found {mp.__version__}"
        )
    mp.mp.dps = WORKING_DIGITS
    outputs = build_outputs()
    output_directory = Path(__file__).resolve().parent / "generated"

    if arguments.check:
        failures = []
        actual_names = (
            {path.name for path in output_directory.iterdir() if path.is_file()}
            if output_directory.is_dir()
            else set()
        )
        expected_names = set(outputs)
        for name, contents in outputs.items():
            path = output_directory / name
            if not path.is_file():
                failures.append(f"missing: {path}")
            elif path.read_bytes() != contents:
                failures.append(f"differs: {path}")
        for name in sorted(actual_names - expected_names):
            failures.append(f"unexpected: {output_directory / name}")
        if failures:
            raise SystemExit("\n".join(failures))
        print(f"verified {len(outputs)} generated files byte-for-byte")
        return

    output_directory.mkdir(parents=True, exist_ok=True)
    for name, contents in outputs.items():
        (output_directory / name).write_bytes(contents)
    print(f"wrote {len(outputs)} files to {output_directory}")


if __name__ == "__main__":
    main()
