#!/usr/bin/env python3
"""Stage 7K.2: Out-of-distribution geometry-coordinate validation.

Fit on Stage 7J.3 constant-area geometry perturbations.
Test on Stage 7J independently varied outlet radius and width.

No nonlinear solver calls. Standard library only.
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


R0 = 0.20
B0 = 0.025
LAMBDA0 = B0 / R0
N0 = 12000.0
NSCALE = 3000.0


def mean(values):
    return sum(values) / len(values)


def rmse(values):
    return math.sqrt(mean([x*x for x in values]))


def solve(a, b):
    n = len(b)
    m = [list(a[i]) + [b[i]] for i in range(n)]

    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(m[i][k]))
        m[k], m[pivot] = m[pivot], m[k]

        if abs(m[k][k]) < 1e-14:
            raise ValueError("Singular regression")

        d = m[k][k]

        for j in range(k, n+1):
            m[k][j] /= d

        for i in range(n):
            if i == k:
                continue

            d = m[i][k]

            for j in range(k, n+1):
                m[i][j] -= d*m[k][j]

    return [m[i][-1] for i in range(n)]


def fit(x, y):
    n = len(x[0])

    gram = [
        [
            sum(row[i]*row[j] for row in x)
            for j in range(n)
        ]
        for i in range(n)
    ]

    rhs = [
        sum(row[i]*v for row, v in zip(x, y))
        for i in range(n)
    ]

    return solve(gram, rhs)


def load(path, geometry_mode):
    rows = []

    with path.open(newline="") as f:
        for record in csv.DictReader(f):
            if record["classification"] != "interior_peak":
                continue

            rpm = float(record["rpm"])
            r = float(record["outlet_radius"])
            b = float(record["outlet_width"])
            area = float(record["outlet_effective_area"])
            mdot = float(record["peak_flow"])

            if geometry_mode == "constant":
                tip = float(record["tip_speed"])
                family = "constant_area"
            else:
                tip = 2*math.pi*rpm*r/60
                family = record["variable"]

            q = mdot/(area*tip)

            rows.append({
                "rpm": rpm,
                "scale": float(record["scale"]),
                "family": family,
                "r": r,
                "b": b,
                "area": area,
                "tip": tip,
                "mdot": mdot,
                "q": q,
                "x": math.log((b/r)/LAMBDA0),
                "z": (rpm-N0)/NSCALE,
            })

    return rows


def add_baselines(training, testing):
    baseline = {}

    for row in training:
        if abs(row["scale"]-1) < 1e-12:
            baseline[row["rpm"]] = row["q"]

    for row in training + testing:
        if row["rpm"] not in baseline:
            raise ValueError("Missing baseline speed")

        row["baseline_q"] = baseline[row["rpm"]]
        row["y"] = math.log(row["q"]/row["baseline_q"])


def alpha(row, coefficients):
    return coefficients[0] + coefficients[1]*row["z"]


def residual(row, coefficients):
    return row["y"] + alpha(row, coefficients)*row["x"]


def metrics(rows, coefficients):
    values = [residual(row, coefficients) for row in rows]

    return {
        "count": len(values),
        "rmse": rmse(values),
        "bias": mean(values),
        "max_abs": max(abs(x) for x in values),
    }


def print_metrics(label, rows, coefficients):
    m = metrics(rows, coefficients)

    print(
        f"{label:32s} "
        f"n={m['count']:3d} "
        f"RMSE={m['rmse']:.9e} "
        f"bias={m['bias']:+.9e} "
        f"max={m['max_abs']:.9e}"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--constant",
        type=Path,
        default=Path(
            "benchmarks/results-constant-area-convergence/"
            "constant_area_sensitivity.csv"
        )
    )

    parser.add_argument(
        "--independent",
        type=Path,
        default=Path(
            "benchmarks/results-geometry-sensitivity/"
            "geometry_sensitivity.csv"
        )
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "benchmarks/results-stage7k2/geometry_transfer.csv"
        )
    )

    args = parser.parse_args()

    training = load(args.constant, "constant")
    testing = load(args.independent, "independent")

    add_baselines(training, testing)

    # Fit exactly the Stage 7K.1 speed-dependent exponent:
    #
    # log(Q/Q0) = -[a+b*z] log(lambda/lambda0)

    nonbaseline = [
        row for row in training
        if abs(row["x"]) > 1e-12
    ]

    coefficients = fit(
        [[row["x"], row["x"]*row["z"]] for row in nonbaseline],
        [-row["y"] for row in nonbaseline]
    )

    print("STAGE 7K.2 GEOMETRY TRANSFER")
    print("="*90)

    print(f"Constant-area training cases: {len(training)}")
    print(f"Independent-geometry test cases: {len(testing)}")

    print(
        "\nFitted exponent: "
        f"alpha={coefficients[0]:+.9f}"
        f"{coefficients[1]:+.9f}*(RPM-12000)/3000"
    )

    print("\nTRAINING GEOMETRY")
    print("-"*90)

    print_metrics(
        "Uncorrected",
        training,
        [0.0, 0.0]
    )

    print_metrics(
        "Speed-dependent correction",
        training,
        coefficients
    )

    print("\nOUT-OF-DISTRIBUTION GEOMETRY")
    print("-"*90)

    print_metrics(
        "Uncorrected",
        testing,
        [0.0, 0.0]
    )

    print_metrics(
        "Transferred correction",
        testing,
        coefficients
    )

    for family in ("outlet_radius", "outlet_width"):
        subset = [
            row for row in testing
            if row["family"] == family
        ]

        print(f"\n{family}")
        print_metrics("Uncorrected", subset, [0.0, 0.0])
        print_metrics("Transferred correction", subset, coefficients)

    print("\nBY SPEED AND GEOMETRY FAMILY")
    print("-"*90)

    groups = defaultdict(list)

    for row in testing:
        groups[(row["rpm"], row["family"])].append(row)

    for (rpm, family), rows in sorted(groups.items()):
        before = metrics(rows, [0.0, 0.0])
        after = metrics(rows, coefficients)

        print(
            f"RPM={rpm:6.0f} {family:14s} "
            f"before={before['rmse']:.9e} "
            f"after={after['rmse']:.9e}"
        )

    # Diagnostic only: fit separate log-area and log-aspect
    # corrections on the independent geometry observations.
    #
    # This is an in-sample fit, NOT transfer validation.
    # It estimates whether an additional geometry direction
    # is needed.

    a0 = 2*math.pi*R0*B0*(1-0.10)

    diagnostic = [
        row for row in testing
        if abs(row["scale"]-1) > 1e-12
    ]

    x = [
        [
            math.log(row["area"]/a0),
            row["x"]
        ]
        for row in diagnostic
    ]

    y = [-row["y"] for row in diagnostic]

    area_coefficient, aspect_coefficient = fit(x, y)

    diagnostic_errors = [
        row["y"]
        + area_coefficient*math.log(row["area"]/a0)
        + aspect_coefficient*row["x"]
        for row in diagnostic
    ]

    print("\nIN-SAMPLE TWO-DIRECTION DIAGNOSTIC")
    print("-"*90)

    print(
        "log-area coefficient:   "
        f"{area_coefficient:+.9f}"
    )

    print(
        "log-aspect coefficient: "
        f"{aspect_coefficient:+.9f}"
    )

    print(
        "In-sample log-RMSE:     "
        f"{rmse(diagnostic_errors):.9e}"
    )

    print(
        "\nThis diagnostic is not a validated "
        "predictive model."
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", newline="") as f:
        fields = [
            "rpm", "family", "scale",
            "radius", "width", "area", "tip_speed",
            "peak_flow", "qstar", "baseline_qstar",
            "log_aspect_ratio", "alpha",
            "uncorrected_log_error",
            "corrected_log_error",
            "corrected_qstar"
        ]

        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()

        for row in testing:
            a = alpha(row, coefficients)

            writer.writerow({
                "rpm": row["rpm"],
                "family": row["family"],
                "scale": row["scale"],
                "radius": row["r"],
                "width": row["b"],
                "area": row["area"],
                "tip_speed": row["tip"],
                "peak_flow": row["mdot"],
                "qstar": row["q"],
                "baseline_qstar": row["baseline_q"],
                "log_aspect_ratio": row["x"],
                "alpha": a,
                "uncorrected_log_error": row["y"],
                "corrected_log_error":
                    residual(row, coefficients),
                "corrected_qstar":
                    row["q"]*math.exp(a*row["x"]),
            })

    print(f"\nOutput: {args.output}")


if __name__ == "__main__":
    main()
