#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    df = pd.read_csv(args.csv)

    output = args.output or args.csv.parent
    output.mkdir(parents=True, exist_ok=True)

    temperatures = sorted(df.initial_T.unique())
    pressures = sorted(df.initial_p.unique())

    for model in sorted(df.model.unique()):
        for loss in sorted(df.loss.unique()):
            fig, axes = plt.subplots(
                1, 2, figsize=(12, 5),
                sharex=True, sharey=True
            )

            for ax, solver in zip(axes, ["full", "fast"]):
                subset = df[
                    (df.model == model) &
                    (df.loss == loss) &
                    (df.solver == solver)
                ]

                grid = np.full(
                    (len(temperatures), len(pressures)),
                    np.nan
                )

                for _, row in subset.iterrows():
                    i = temperatures.index(row.initial_T)
                    j = pressures.index(row.initial_p)

                    if row.recovered:
                        grid[i, j] = row.iterations
                    elif row.converged:
                        grid[i, j] = -1
                    else:
                        grid[i, j] = -2

                masked = np.ma.masked_where(grid < 0, grid)

                image = ax.imshow(
                    masked,
                    origin="lower",
                    aspect="auto",
                    interpolation="nearest",
                    cmap="viridis",
                    vmin=0,
                    vmax=40
                )

                for i in range(len(temperatures)):
                    for j in range(len(pressures)):
                        if grid[i, j] == -1:
                            ax.text(
                                j, i, "O",
                                ha="center", va="center",
                                fontsize=8
                            )
                        elif grid[i, j] == -2:
                            ax.text(
                                j, i, "X",
                                ha="center", va="center",
                                fontsize=8
                            )

                ax.set_title(f"{solver} solver")
                ax.set_xticks(range(len(pressures)))
                ax.set_xticklabels(
                    [f"{p/1000:g}" for p in pressures],
                    rotation=45
                )
                ax.set_yticks(range(len(temperatures)))
                ax.set_yticklabels(
                    [f"{t:g}" for t in temperatures]
                )
                ax.set_xlabel("Initial pressure (kPa)")
                ax.set_ylabel("Initial temperature (K)")

            fig.colorbar(
                image, ax=axes,
                label="Newton iterations",
                shrink=0.8
            )

            fig.suptitle(
                f"{model}: entropy generation = {loss:g} J/(kg K)\n"
                "X = failed solve, O = converged to another root"
            )

            filename = (
                f"basin_{model}_loss{loss:g}.png"
            )

            fig.savefig(
                output / filename,
                dpi=180,
                bbox_inches="tight"
            )

            plt.close(fig)

            print(output / filename)


if __name__ == "__main__":
    main()
