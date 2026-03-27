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

from utils import (
    METHODS,
    get_rmsd_dataframe,
    plateau_from_cage_size,
    get_diffusion_coefficient,
)

sns.set(font_scale=2)


def main(folder: str) -> None:
    files = glob.glob(f"{folder}/*/*/*.csv")

    df_list = []
    for f in files:
        ddf = pd.read_csv(f)
        ddf["rate"] = float(f.split("/")[-3])
        df_list.append(ddf)

    df = pd.concat(df_list).reset_index(drop=True)

    id_columns = ["time", "rate"]

    df_rmsd = get_rmsd_dataframe(df, id_columns)

    cage_size = 0.1
    plateau = plateau_from_cage_size(cage_size)

    # Largest rate, look at all definitions
    rate = df_rmsd["rate"].max()
    rate = 10000
    ax = sns.lineplot(
        data=df_rmsd[df_rmsd["rate"] == rate],
        x="time",
        y="RMSD",
        hue="Definition",
        linewidth=4,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(y=np.sqrt(plateau), color="black", linestyle="dashed")
    plt.title(f"Escape - rate = {rate}")
    plt.show()

    # All rates, unbounded
    ax = sns.lineplot(
        data=df_rmsd.query("Definition == 'Unbounded'"),
        x="time",
        y="RMSD",
        hue="rate",
        linewidth=4,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(y=np.sqrt(plateau), color="black", linestyle="dashed")
    plt.title("Escape")
    plt.show()

    df_D = (
        df_rmsd.groupby(["rate", "Definition"])
        .apply(get_diffusion_coefficient)
        .rename("D")
        .reset_index()
    )

    x = np.logspace(1, 6, 100)
    y = 0.75e-2 / x

    ax = sns.scatterplot(
        data=df_D[df_D["Definition"].isin(["Integral", "Unbounded"])],
        x="rate",
        y="D",
        hue="Definition",
    )
    ax.plot(x, y, linestyle="dashed", color="grey")
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.show()


if __name__ == "__main__":
    folder = "../production/escape/"
    main(folder)
