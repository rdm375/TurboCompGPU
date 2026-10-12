#!/usr/bin/env python3
"""Stage 7K.1: Geometry-informed compressor peak coordinates.

Analyze constant-area geometry perturbations from Stage 7J.3.

Models:
    Q* = m_peak / (A_eff * U2)

    Q_corrected = Q* * (lambda/lambda0)**alpha

    alpha = a                           constant model
    alpha = a + b*(rpm-12000)/3000      speed-dependent model

No nonlinear operating-point calculations are performed.
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


BASELINE_RADIUS = 0.2
BASELINE_WIDTH = 0.025
BASELINE_ASPECT = BASELINE_WIDTH / BASELINE_RADIUS


def mean(values):
    return sum(values) / len(values)


def rmse(values):
    return math.sqrt(mean([v * v for v in values]))


def solve_linear(matrix, rhs):
    """Small dense linear system with partial pivoting."""
    n = len(rhs)
    a = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]

    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(a[i][k]))

        if abs(a[pivot][k]) < 1e-14:
            raise ValueError("Singular regression system")

        a[k], a[pivot] = a[pivot], a[k]

        scale = a[k][k]

        for j in range(k, n + 1):
            a[k][j] /= scale

        for i in range(n):
            if i == k:
                continue

            factor = a[i][k]

            for j in range(k, n + 1):
                a[i][j] -= factor * a[k][j]

    return [a[i][n] for i in range(n)]


def least_squares(features, targets):
    n = len(features[0])

    gram = [
        [
            sum(row[i] * row[j] for row in features)
            for j in range(n)
        ]
        for i in range(n)
    ]

    rhs = [
        sum(row[i] * y for row, y in zip(features, targets))
        for i in range(n)
    ]

    return solve_linear(gram, rhs)


def speed_coordinate(rpm):
    return (rpm - 12000.0) / 3000.0


def load_data(path):
    data = []

    with path.open(newline="") as f:
        reader = csv.DictReader(f)

        required = {
            "rpm",
            "scale",
            "outlet_radius",
            "outlet_width",
            "outlet_effective_area",
            "tip_speed",
            "classification",
            "peak_flow",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(f"Missing CSV columns: {sorted(missing)}")

        for row in reader:
            rpm = float(row["rpm"])
            scale = float(row["scale"])
            radius = float(row["outlet_radius"])
            width = float(row["outlet_width"])
            area = float(row["outlet_effective_area"])
            tip_speed = float(row["tip_speed"])

            classification = row["classification"]

            if classification != "interior_peak":
                continue

            peak_flow = float(row["peak_flow"])

            if not all(
                math.isfinite(v) and v > 0
                for v in (radius, width, area, tip_speed, peak_flow)
            ):
                raise ValueError("Invalid positive geometry or peak value")

            aspect = width / radius

            qstar = peak_flow / (area * tip_speed)

            data.append({
                "rpm": rpm,
                "scale": scale,
                "radius": radius,
                "width": width,
                "area": area,
                "tip_speed": tip_speed,
                "peak_flow": peak_flow,
                "aspect": aspect,
                "qstar": qstar,
                "log_aspect": math.log(aspect / BASELINE_ASPECT),
            })

    return data


def attach_baselines(data):
    baselines = {}

    for row in data:
        if abs(row["scale"] - 1.0) < 1e-12:
            baselines[row["rpm"]] = row["qstar"]

    for row in data:
        if row["rpm"] not in baselines:
            raise ValueError(
                f"Missing baseline for RPM={row['rpm']}"
            )

        row["baseline_qstar"] = baselines[row["rpm"]]

        row["log_response"] = math.log(
            row["qstar"] / row["baseline_qstar"]
        )


def design_row(row, model):
    x = row["log_aspect"]

    if model == "constant":
        return [x]

    if model == "speed":
        z = speed_coordinate(row["rpm"])
        return [x, x * z]

    raise ValueError(model)


def fit(data, model):
    training = [
        row for row in data
        if abs(row["log_aspect"]) > 1e-12
    ]

    if not training:
        raise ValueError("No nonbaseline training observations")

    features = [design_row(row, model) for row in training]

    # log(Q*/Q*_baseline) = -alpha * log(lambda/lambda0)
    targets = [-row["log_response"] for row in training]

    return least_squares(features, targets)


def exponent(row, coefficients, model):
    if model == "constant":
        return coefficients[0]

    return (
        coefficients[0]
        + coefficients[1] * speed_coordinate(row["rpm"])
    )


def corrected_log_response(row, coefficients, model):
    alpha = exponent(row, coefficients, model)

    return row["log_response"] + alpha * row["log_aspect"]


def evaluate(data, coefficients, model):
    errors = [
        corrected_log_response(row, coefficients, model)
        for row in data
        if abs(row["log_aspect"]) > 1e-12
    ]

    return rmse(errors) if errors else math.nan


def baseline_error(data):
    errors = [
        row["log_response"]
        for row in data
        if abs(row["log_aspect"]) > 1e-12
    ]

    return rmse(errors) if errors else math.nan


def print_model_comparison(data):
    print("\nFULL-DATA FIT")
    print("=" * 78)

    constant = fit(data, "constant")
    speed = fit(data, "speed")

    print(f"Uncorrected log-RMSE: {baseline_error(data):.9e}")
    print(f"Constant-alpha log-RMSE: {evaluate(data, constant, 'constant'):.9e}")
    print(f"Speed-alpha log-RMSE: {evaluate(data, speed, 'speed'):.9e}")

    print(f"\nConstant alpha: {constant[0]:+.9f}")

    print(
        "Speed-dependent alpha: "
        f"{speed[0]:+.9f} "
        f"{speed[1]:+.9f}*(RPM-12000)/3000"
    )

    return constant, speed


def leave_one_speed_out(data):
    print("\nLEAVE-ONE-SPEED-OUT VALIDATION")
    print("=" * 78)

    print(
        f"{'RPM':>8}"
        f"{'uncorrected':>18}"
        f"{'constant':>18}"
        f"{'speed-dependent':>20}"
    )

    pooled = {
        "uncorrected": [],
        "constant": [],
        "speed": [],
    }

    for rpm in sorted({row["rpm"] for row in data}):
        training = [row for row in data if row["rpm"] != rpm]
        testing = [
            row for row in data
            if row["rpm"] == rpm
            and abs(row["log_aspect"]) > 1e-12
        ]

        c = fit(training, "constant")
        s = fit(training, "speed")

        e0 = [row["log_response"] for row in testing]

        ec = [
            corrected_log_response(row, c, "constant")
            for row in testing
        ]

        es = [
            corrected_log_response(row, s, "speed")
            for row in testing
        ]

        pooled["uncorrected"].extend(e0)
        pooled["constant"].extend(ec)
        pooled["speed"].extend(es)

        print(
            f"{rpm:8.0f}"
            f"{rmse(e0):18.9e}"
            f"{rmse(ec):18.9e}"
            f"{rmse(es):20.9e}"
        )

    print("\nPooled held-out log-RMSE:")

    for key, errors in pooled.items():
        print(f"  {key:18s} {rmse(errors):.9e}")


def held_out_geometry(data):
    print("\nHELD-OUT GEOMETRY VALIDATION")
    print("=" * 78)

    # Fit only near baseline.
    # Test on larger geometry changes.
    #
    # The baseline is always available at each speed.
    training = [
        row for row in data
        if abs(row["scale"] - 1.0) <= 0.01 + 1e-12
    ]

    testing = [
        row for row in data
        if abs(row["scale"] - 1.0) >= 0.02 - 1e-12
    ]

    c = fit(training, "constant")
    s = fit(training, "speed")

    print(f"Training cases: {len(training)}")
    print(f"Testing cases:  {len(testing)}")

    print(f"Uncorrected:      {baseline_error(testing):.9e}")
    print(f"Constant alpha:   {evaluate(testing, c, 'constant'):.9e}")
    print(f"Speed alpha:      {evaluate(testing, s, 'speed'):.9e}")

    print(f"\nConstant alpha: {c[0]:+.9f}")

    print(
        "Speed alpha: "
        f"{s[0]:+.9f} "
        f"{s[1]:+.9f}*(RPM-12000)/3000"
    )


def write_output(path, data, constant, speed):
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "rpm",
        "scale",
        "radius",
        "width",
        "area",
        "tip_speed",
        "peak_flow",
        "aspect",
        "qstar",
        "baseline_qstar",
        "log_response",
        "qstar_constant_corrected",
        "qstar_speed_corrected",
        "alpha_constant",
        "alpha_speed",
    ]

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()

        for row in data:
            ac = exponent(row, constant, "constant")
            a_speed = exponent(row, speed, "speed")

            record = {
                key: row[key]
                for key in fields
                if key in row
            }

            record.update({
                "qstar_constant_corrected":
                    row["qstar"] * math.exp(
                        ac * row["log_aspect"]
                    ),
                "qstar_speed_corrected":
                    row["qstar"] * math.exp(
                        a_speed * row["log_aspect"]
                    ),
                "alpha_constant": ac,
                "alpha_speed": a_speed,
            })

            writer.writerow(record)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input",
        type=Path,
        help="Stage 7J.3 convergence CSV"
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "benchmarks/results-stage7k1/peak_coordinates.csv"
        )
    )

    args = parser.parse_args()

    data = load_data(args.input)

    attach_baselines(data)

    print("STAGE 7K.1 PEAK COORDINATE ANALYSIS")
    print("=" * 78)

    print(f"Interior peak cases: {len(data)}")
    print(f"Rotational speeds: {len(set(r['rpm'] for r in data))}")

    grouped = defaultdict(list)

    for row in data:
        grouped[row["rpm"]].append(row)

    print("\nRAW PEAK COORDINATE SPREAD")
    print("=" * 78)

    for rpm, rows in sorted(grouped.items()):
        q = [row["qstar"] for row in rows]

        spread = (max(q) - min(q)) / mean(q)

        print(
            f"RPM={rpm:6.0f} "
            f"cases={len(rows):2d} "
            f"relative_Qstar_spread={spread:.9e}"
        )

    constant, speed = print_model_comparison(data)

    leave_one_speed_out(data)

    held_out_geometry(data)

    write_output(args.output, data, constant, speed)

    print(f"\nOutput: {args.output}")


if __name__ == "__main__":
    main()
