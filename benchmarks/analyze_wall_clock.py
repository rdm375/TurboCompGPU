#!/usr/bin/env python3

import sys
from pathlib import Path

import numpy as np
import pandas as pd

path = Path(
    sys.argv[1]
    if len(sys.argv) > 1
    else "benchmarks/results-wall-clock/timings.csv"
)

df = pd.read_csv(path)

required = {
    "model", "loss", "mode",
    "repetition", "elapsed_ns"
}

missing = required - set(df.columns)

if missing:
    raise ValueError(f"Missing CSV columns: {sorted(missing)}")

if df.empty:
    raise ValueError("Benchmark CSV is empty")

if not np.isfinite(df.elapsed_ns).all():
    raise ValueError("Non-finite timing values")

if (df.elapsed_ns <= 0).any():
    raise ValueError("Nonpositive timing values")

keys = ["model", "loss"]

summary = (
    df.groupby(keys + ["mode"])
    .elapsed_ns
    .agg(
        samples="count",
        median_ns="median",
        mean_ns="mean",
        minimum_ns="min",
        p05_ns=lambda x: x.quantile(0.05),
        p95_ns=lambda x: x.quantile(0.95),
    )
    .reset_index()
)

for column in [
    "median_ns", "mean_ns", "minimum_ns",
    "p05_ns", "p95_ns"
]:
    summary[column.replace("_ns", "_us")] = (
        summary[column] / 1000.0
    )

median = summary.pivot(
    index=keys,
    columns="mode",
    values="median_ns"
)

expected_modes = {"full", "fast", "fast_validated"}

if set(median.columns) != expected_modes:
    raise ValueError(
        f"Unexpected benchmark modes: {list(median.columns)}"
    )

speedup = pd.DataFrame(index=median.index)

speedup["full_us"] = median["full"] / 1000.0
speedup["fast_us"] = median["fast"] / 1000.0
speedup["validated_us"] = median["fast_validated"] / 1000.0

speedup["fast_speedup"] = (
    median["full"] / median["fast"]
)

speedup["validated_speedup"] = (
    median["full"] / median["fast_validated"]
)

speedup["validation_overhead_us"] = (
    median["fast_validated"] - median["fast"]
) / 1000.0

print("\nTIMING DISTRIBUTIONS (microseconds)\n")

print(
    summary[
        [
            "model", "loss", "mode", "samples",
            "minimum_us", "p05_us",
            "median_us", "p95_us"
        ]
    ].to_string(index=False, float_format=lambda x: f"{x:.4f}")
)

print("\nMEDIAN-BASED SPEEDUPS\n")

print(
    speedup.to_string(
        float_format=lambda x: f"{x:.4f}"
    )
)

# Additional paired ratios compare the same repetition index.
# These are useful alongside ratios of marginal medians.
paired = df.pivot(
    index=keys + ["repetition"],
    columns="mode",
    values="elapsed_ns"
)

if paired.isna().any().any():
    raise ValueError("Incomplete paired timing observations")

paired["fast_ratio"] = paired["full"] / paired["fast"]
paired["validated_ratio"] = (
    paired["full"] / paired["fast_validated"]
)

paired_summary = (
    paired.groupby(level=keys)[
        ["fast_ratio", "validated_ratio"]
    ]
    .median()
)

print("\nMEDIAN PAIRED SPEEDUP RATIOS\n")

print(
    paired_summary.to_string(
        float_format=lambda x: f"{x:.4f}"
    )
)

output_dir = path.parent

summary.to_csv(
    output_dir / "timing_summary.csv",
    index=False
)

speedup.to_csv(
    output_dir / "speedup_summary.csv"
)

paired_summary.to_csv(
    output_dir / "paired_speedup_summary.csv"
)

print("\nPASS: Timing analysis completed")
