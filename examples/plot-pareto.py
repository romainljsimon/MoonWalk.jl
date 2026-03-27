# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pandas>=3.0.1",
#     "seaborn>=0.13.2",
# ]
# ///

import glob
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set(font_scale=2)

methods = ["ExactRotation", "Integral", "Unbounded"]


def main(folder: str) -> None:
    files = glob.glob(f"{folder}/*/*/*.csv")

    df_list = []
    for f in files:
        ddf = pd.read_csv(f)
        ddf["alpha"] = f.split("/")[1]
        df_list.append(ddf)

    df = pd.concat(df_list).reset_index(drop=True)

    # assert len(df) == len(files) * 56

    for method in methods:
        df[f"{method}"] = (
            df[f"{method}_x"] ** 2 + df[f"{method}_y"] ** 2 + df[f"{method}_z"] ** 2
        )

    df_rmsd = (
        df[["time", "alpha"] + methods]
        .melt(id_vars=["time", "alpha"], value_vars=methods, var_name="Definition")
        .groupby(["time", "alpha", "Definition"])
        .apply(lambda x: np.sqrt(np.mean(x)))
        .rename("RMSD")
        .reset_index()
    )

    # Diffusive with a dummy coefficient
    x = np.logspace(1, 5, 1000)
    y = np.sqrt(0.05 * x)

    for definition in ["Unbounded", "Integral"]:
        ax = sns.lineplot(
            data=df_rmsd.query(f"Definition == '{definition}'"),
            x="time",
            y="RMSD",
            hue="alpha",
            linewidth=4,
        )
        ax.plot(x, y, linestyle="dashed", color="grey")
        ax.set_xscale("log")
        ax.set_yscale("log")
        plt.title(f"CTRW - {definition}")
        plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} folder")
    main(sys.argv[1])
