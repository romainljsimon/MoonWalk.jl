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
import scipy

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
        ddf["alpha"] = float(f.split("/")[-3])
        df_list.append(ddf)

    df = pd.concat(df_list).reset_index(drop=True)

    alphas = set(df["alpha"])

    times = sorted([1e7, 1e6, 1e5, 1e4])
    ddf = df[(df["time"].isin(times))][
        ["time", "alpha", "Unbounded_x", "Unbounded_y", "Unbounded_z"]
    ].melt(id_vars=["time", "alpha"])[["time", "alpha", "value"]]

    print(ddf.query("alpha == 0.7").value_counts("time"))
    print(ddf.query("alpha == 1.2").value_counts("time"))

    def compute_distrib(df, n_bins=40):
        min_v, max_v = np.percentile(df["value"], [0.5, 99.5])
        df = df[(df["value"] > min_v) & (df["value"] < max_v)]

        density, bin_edges = np.histogram(df["value"], bins=30, density=True)
        bin_centers = [
            (bin_edges[i + 1] + bin_edges[i]) / 2 for i in range(len(bin_edges) - 1)
        ]

        return pd.DataFrame({"angle": bin_centers, "density": density})

    df_distrib = (
        ddf.groupby(["time", "alpha"])
        .apply(compute_distrib)
        .reset_index()
        .drop(columns="level_2")
    )

    for time in times:
        for alpha in alphas:
            filename = f"distrib_{alpha}_{time:.0e}.csv"
            df_distrib.query(f"time == {time}").query(
                f"alpha == {alpha}"
            ).reset_index()[["angle", "density"]].to_csv(filename, sep=" ", index=False)
            add_pound_key(filename)

    def plot_distrib(data, **kwargs):
        alpha = data.iloc[0]["alpha"]
        time = data.iloc[0]["time"]

        df_local = df_distrib.query(f"alpha == {alpha}").query(f"time == {time}")

        x = np.linspace(df_local["angle"].min(), df_local["angle"].max(), num=1000)
        y = scipy.stats.norm.pdf(x, loc=0, scale=np.std(data["value"]))

        ax = sns.lineplot(data=df_local, x="angle", y="density", linewidth=4)
        ax.plot(x, y, linestyle="dashed", color="grey", linewidth=4)

    def reformat_label(label):
        a, b = label.split("|")
        alpha = float(a.split("=")[1])
        time = float(b.split("=")[1])
        return f"$\\alpha$ = {alpha}, time = {time:2.2g}"

    g = sns.FacetGrid(ddf, col="time", row="alpha", sharex=False)  # , sharey=False)
    # Loop through all axes in the FacetGrid
    for ax in g.axes.flat:
        ax.title.set_text(reformat_label(ax.title.get_text()))
    g.map_dataframe(plot_distrib)
    g.set(yscale="log")
    plt.tight_layout()
    plt.show()

    id_columns = ["time", "alpha"]

    df_msd = get_msd_dataframe(df, id_columns)

    cage_size = 0.2
    plateau = plateau_from_cage_size(cage_size)

    df_msd["$\\alpha$"] = df_msd["alpha"].apply(str)

    # Save in wide format
    for alpha in set(df_msd["alpha"]):
        df_msd.query(f"alpha == {alpha}").pivot(
            index="time", values="MSD", columns="Definition"
        ).reset_index().to_csv(f"alpha_{alpha}_msd.csv", sep=" ", index=False)
        add_pound_key(f"alpha_{alpha}_msd.csv")

    # one alpha at a time, look at all definitions
    x = np.logspace(4, 7, 1000)
    y = 0.5e-5 * x
    for alpha in set(df_msd["alpha"]):
        ax = sns.lineplot(
            data=df_msd[df_msd["alpha"] == alpha],
            x="time",
            y="MSD",
            hue="Definition",
            linewidth=4,
        )
        ax.plot(x, y, linestyle="dashed", color="grey")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(y=plateau, color="black", linestyle="dashed")
        ax.axhline(y=(np.pi**2 + 6) / 3, color="grey", linestyle="dashed")
        plt.title(f"Pareto - $\\alpha$ = {alpha}")
        plt.show()

        D = get_diffusion_coefficient(
            df_msd[(df_msd["alpha"] == alpha) & (df_msd["Definition"] == "Integral")],
            title=f"Pareto - $\\alpha$ = {alpha} - Integral method",
        )


if __name__ == "__main__":
    folder = "../production/pareto/"
    main(folder)
