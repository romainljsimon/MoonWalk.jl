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

    assert len(df) == len(files) * 45

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

    # Theoretical value
    # Jump at times between [0, T], of an amplitude between [0, A]
    # RMSD is sqrt(A**2 * t / T**2)
    T = 1
    A = 0.1 / 2
    x = np.logspace(1, 4, 1000)
    y = np.sqrt(2 * A**2 * x / T)

    ax = sns.lineplot(data=df_rmsd, x="time", y="RMSD", hue="Definition", linewidth=4)
    ax.plot(x, y, linestyle="dashed", color="grey")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(y=np.sqrt((np.pi**2 + 6) / 3), color="grey", linestyle="dashed")
    plt.title("Brownian diffusion")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} folder")
    main(sys.argv[1])
