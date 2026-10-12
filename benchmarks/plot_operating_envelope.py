#!/usr/bin/env python3
"""M1.5 operating-envelope visualization and paired diagnostics."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("benchmarks/results-operating-envelope")
CSV = ROOT / "envelope.csv"

df = pd.read_csv(CSV)

required = {
    "model", "method", "mass_flow", "rotational_speed_rpm",
    "converged", "validated", "outlet_temperature",
    "outlet_pressure", "pressure_ratio", "specific_work",
    "shaft_power", "iterations", "evaluations",
}

missing = required - set(df.columns)

if missing:
    raise RuntimeError(f"Missing CSV columns: {sorted(missing)}")

keys = ["model", "mass_flow", "rotational_speed_rpm"]

a = (
    df[df["method"] == "independent"]
    .set_index(keys)
    .sort_index()
)

b = (
    df[df["method"] == "continuation"]
    .set_index(keys)
    .sort_index()
)

if not a.index.equals(b.index):
    raise RuntimeError("Independent/continuation cases do not match")

if len(a) != 70 or len(b) != 70:
    raise RuntimeError("Unexpected operating-envelope dimensions")

if not ((df["converged"] == 1) & (df["validated"] == 1)).all():
    raise RuntimeError("Some operating points failed convergence/validation")

# Regression classification thresholds, not Newton tolerances.
tolerances = {
    "outlet_temperature": 1.0e-6,
    "outlet_pressure": 1.0e-3,
    "specific_work": 1.0e-5,
}

report = []

for model in sorted(df["model"].unique()):
    left = a.loc[model]
    right = b.loc[model]

    report.append(f"Model: {model}")
    report.append(f"Operating conditions: {len(left)}")

    for quantity, tolerance in tolerances.items():
        delta = np.abs(
            left[quantity].to_numpy() -
            right[quantity].to_numpy()
        )

        maximum = float(delta.max())

        report.append(
            f"Maximum |delta {quantity}| = {maximum:.12e}"
        )

        if maximum > tolerance:
            raise RuntimeError(
                f"{model}: paired disagreement in {quantity}"
            )

    independent_iterations = left["iterations"].sum()
    continuation_iterations = right["iterations"].sum()

    independent_evaluations = left["evaluations"].sum()
    continuation_evaluations = right["evaluations"].sum()

    report.append(
        "Iteration reduction: "
        f"{100 * (1 - continuation_iterations / independent_iterations):.2f}%"
    )

    report.append(
        "Evaluation reduction: "
        f"{100 * (1 - continuation_evaluations / independent_evaluations):.2f}%"
    )

    report.append("")

    independent = df[
        (df["model"] == model) &
        (df["method"] == "independent")
    ].copy()

    quantities = [
        ("pressure_ratio", "Pressure ratio", "pressure_ratio"),
        ("specific_work", "Specific work (J/kg)", "specific_work"),
        ("shaft_power", "Shaft power (W)", "shaft_power"),
    ]

    for column, ylabel, filename in quantities:
        fig, ax = plt.subplots(figsize=(8, 5))

        for speed, group in independent.groupby(
            "rotational_speed_rpm"
        ):
            group = group.sort_values("mass_flow")

            ax.plot(
                group["mass_flow"],
                group[column],
                marker="o",
                label=f"{speed:g} rpm",
            )

        ax.set_xlabel("Mass flow (kg/s)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{model}: {ylabel}")
        ax.grid(True, alpha=0.3)
        ax.legend(title="Rotational speed")
        fig.tight_layout()

        output = ROOT / f"{model}_{filename}.png"
        fig.savefig(output, dpi=160)
        plt.close(fig)

        print("Created:", output)

    fig, ax = plt.subplots(figsize=(8, 5))

    for method, group in df[df["model"] == model].groupby("method"):
        means = (
            group.groupby("mass_flow")["iterations"]
            .mean()
            .sort_index()
        )

        ax.plot(
            means.index,
            means.values,
            marker="o",
            label=method,
        )

    ax.set_xlabel("Mass flow (kg/s)")
    ax.set_ylabel("Mean Newton iterations")
    ax.set_title(f"{model}: Initialization comparison")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    output = ROOT / f"{model}_iterations.png"
    fig.savefig(output, dpi=160)
    plt.close(fig)

    print("Created:", output)

report.append("PASS: All paired operating-point comparisons")

summary = ROOT / "validation_summary.txt"
summary.write_text("\n".join(report) + "\n", encoding="utf-8")

print()
print(summary.read_text())
