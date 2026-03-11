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
    files = glob.glob(f"{folder}/**/*.csv")

    df = pd.concat([pd.read_csv(f) for f in files]).reset_index(drop=True)

    # assert len(df) == len(files) * 45

    for method in methods:
        df[f"{method}"] = (
            df[f"{method}_x"] ** 2 + df[f"{method}_y"] ** 2 + df[f"{method}_z"] ** 2
        )

    df_rmsd = (
        df[["time"] + methods]
        .melt(id_vars="time", value_vars=methods, var_name="Definition")
        .groupby(["time", "Definition"])
        .apply(lambda x: np.sqrt(np.mean(x)))
        .rename("RMSD")
        .reset_index()
    )

    H = 1
    b = (-(3 * H**2 - 6) * np.sin(H) - 6 * H * np.cos(H) + H**3) / (3 * (H - np.sin(H)))

    ax = sns.lineplot(data=df_rmsd, x="time", y="RMSD", hue="Definition", linewidth=4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(y=np.sqrt(b), color="grey", linestyle="dashed")
    ax.axhline(y=H, color="black", linestyle="dashed")
    plt.title("Cage")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} folder")
    main(sys.argv[1])
