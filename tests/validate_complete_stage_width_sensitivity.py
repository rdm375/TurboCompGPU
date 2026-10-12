#!/usr/bin/env python3

"""Regression validation for Stage 8B.5c.2."""

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: validate_complete_stage_width_sensitivity.py CSV"
        )

    path = Path(sys.argv[1])

    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 60, (
        f"Expected 60 configurations, found {len(rows)}"
    )

    groups = defaultdict(list)
    speed_groups = defaultdict(list)

    for row in rows:
        model = row["model"]
        rpm = float(row["rpm"])
        width = float(row["width_ratio"])

        assert model in {"IdealGas", "PengRobinson"}
        assert row["recovery_converged"] == "1"
        assert row["capacity_converged"] == "1"
        assert row["capacity_reason"] == "Converged"

        recovery = float(row["recovery_mass_kg_s"])
        capacity = float(row["capacity_mass_kg_s"])
        ratio = float(row["recovery_to_capacity_ratio"])

        assert all(
            math.isfinite(value)
            for value in (recovery, capacity, ratio)
        )

        assert 0.0 < recovery < capacity
        assert math.isclose(
            ratio,
            recovery / capacity,
            rel_tol=1e-10,
            abs_tol=1e-12,
        )

        assert (
            abs(float(row["recovery_residual_pa"]))
            <= 1e-3
        )

        assert (
            abs(float(row["capacity_residual_kg_s"]))
            <= 1e-8
        )

        assert (
            abs(float(row["sonic_mach"]) - 1.0)
            <= 1e-8
        )

        groups[(model, rpm)].append(
            (width, recovery, capacity, ratio)
        )

        speed_groups[(model, width)].append(
            (rpm, ratio)
        )

    assert len(groups) == 10
    assert len(speed_groups) == 12

    for key, values in groups.items():
        values.sort()

        assert len(values) == 6, key

        for previous, current in zip(
            values[:-1], values[1:]
        ):
            assert current[1] > previous[1], (
                "Recovery not increasing with width", key
            )

            assert current[2] > previous[2], (
                "Capacity not increasing with width", key
            )

            assert current[3] < previous[3], (
                "Chi not decreasing with width", key
            )

    for key, values in speed_groups.items():
        values.sort()

        assert len(values) == 5, key

        for previous, current in zip(
            values[:-1], values[1:]
        ):
            assert current[1] > previous[1], (
                "Chi not increasing with speed", key
            )

    print("PASS: 60 converged configurations")
    print("PASS: recovery and capacity residuals")
    print("PASS: sonic critical states")
    print("PASS: monotonic width responses")
    print("PASS: monotonic speed responses")


if __name__ == "__main__":
    main()
