# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pandas>=3.0.1",
#     "scipy>=1.17.1",
#     "seaborn>=0.13.2",
# ]
# ///

import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from utils import (
    get_msd_dataframe,
    plateau_from_cage_size,
    get_diffusion_coefficient,
    add_pound_key,
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

    df_msd = get_msd_dataframe(df, id_columns)

    df_msd.query("Definition == 'Unbounded'").pivot(
        index="time", values="MSD", columns="rate"
    ).reset_index().to_csv("escape_unbounded_msd.csv", sep=" ", index=False)
    add_pound_key("escape_unbounded_msd.csv")

    cage_size = 0.2
    plateau = plateau_from_cage_size(cage_size)

    # Largest rate, look at all definitions
    rate = df_msd["rate"].max()
    rate = 10000
    ax = sns.lineplot(
        data=df_msd[df_msd["rate"] == rate],
        x="time",
        y="MSD",
        hue="Definition",
        linewidth=4,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(y=plateau, color="black", linestyle="dashed")
    plt.title(f"Escape - rate = {rate}")
    plt.show()

    # All rates, unbounded
    ax = sns.lineplot(
        data=df_msd.query("Definition == 'Unbounded'"),
        x="time",
        y="MSD",
        hue="rate",
        linewidth=4,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(y=plateau, color="black", linestyle="dashed")
    plt.title("Escape")
    plt.show()

    df_D = (
        df_msd.groupby(["rate", "Definition"])
        .apply(get_diffusion_coefficient)
        .rename("D")
        .reset_index()
    )

    df_D.pivot(index="rate", values="D", columns="Definition").reset_index()[
        ["rate", "Integral", "Unbounded"]
    ].to_csv("escape_D.csv", sep=" ", index=False)
    add_pound_key("escape_D.csv")

    x = np.logspace(2, 7, 100)
    y = 0.25e-1 / x

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
