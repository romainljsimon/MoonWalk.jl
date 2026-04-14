# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pandas>=3.0.1",
#     "seaborn>=0.13.2",
# ]
# ///

import glob

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
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

    alpha = 1.2
    # times = [1000.0, 43287.612810830615, 10000000.0]
    times = [1e8, 1e7, 1e6, 1e5]
    for time in times:
        ddf = df[(df["alpha"] == alpha) & (df["time"] == time)][
            ["Unbounded_x", "Unbounded_y", "Unbounded_z"]
        ].melt()

        x = np.linspace(ddf["value"].min(), ddf["value"].max(), num=1000)
        y = scipy.stats.norm.pdf(x, loc=0, scale=np.std(ddf["value"]))

        # ax = sns.histplot(
        #    data=ddf,
        #    x="value",
        #    stat="density",
        #    bins=40,
        #    element="step",
        #    fill=False,
        #    linewidth=4,
        # )

        density, bin_edges, _ = ax.hist(ddf["value"], bins=30, density=True)
        bin_centers = [
            (bin_edges[i + 1] + bin_edges[i]) / 2 for i in range(len(bin_edges) - 1)
        ]
        plt.close()

        df2 = pd.DataFrame({"angle": bin_centers, "density": density})

        ax = sns.lineplot(data=df2, x="angle", y="density", linewidth=4)
        ax.plot(x, y, linestyle="dashed", color="grey", linewidth=4)
        ax.set_yscale("log")
        plt.title(f"$\\alpha$ = {alpha}, time = {time:2.2g}")
        plt.xlabel("Angle")
        plt.show()
        ax.set_yscale("log")
        plt.show()

    times = sorted([1e8, 1e7, 1e6, 1e5])
    ddf = df[(df["time"].isin(times))][
        ["time", "alpha", "Unbounded_x", "Unbounded_y", "Unbounded_z"]
    ].melt(id_vars=["time", "alpha"])[["time", "alpha", "value"]]

    def plot_distrib(data, **kwargs):
        x = np.linspace(data["value"].min(), data["value"].max(), num=1000)
        y = scipy.stats.norm.pdf(x, loc=0, scale=np.std(data["value"]))

        density, bin_edges = np.histogram(data["value"], bins=30, density=True)
        bin_centers = [
            (bin_edges[i + 1] + bin_edges[i]) / 2 for i in range(len(bin_edges) - 1)
        ]

        df2 = pd.DataFrame({"angle": bin_centers, "density": density})

        ax = sns.lineplot(data=df2, x="angle", y="density", linewidth=4)
        ax.plot(x, y, linestyle="dashed", color="grey", linewidth=4)
        # ax.set_yscale("log")
        # plt.title(f"$\\alpha$ = {alpha}, time = {time:2.2g}")
        # plt.xlabel("Angle")
        # ax.set_yscale("log")
        # plt.show()

    def reformat_label(label):
        a, b = label.split("|")
        alpha = float(a.split("=")[1])
        time = float(b.split("=")[1])
        return f"$\\alpha$ = {alpha}, time = {time:2.2g}"

    g = sns.FacetGrid(ddf, col="time", row="alpha", sharex=False)
    # Loop through all axes in the FacetGrid
    for ax in g.axes.flat:
        ax.title.set_text(reformat_label(ax.title.get_text()))
    # for t, time in zip(g.axes[0], times):
    #    t.title.set_text(f"time = {time:2.2g}")
    # for t, time in zip(g.axes[1], times):
    #    t.title.set_text(f"time = {time:2.2g}")
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
            plot=True,
            title=f"Pareto - $\\alpha$ = {alpha} - Integral method",
        )


def plot_pareto_cumulative(alpha: float, tau: float) -> None:
    tau = 1000
    x = np.linspace(tau, 20 * tau, num=1000)

    df_list = []
    for alpha in [0.5, 0.8, 1.5]:
        prob = alpha * tau**alpha / (x ** (1 + alpha))
        cumulative = np.array([sum(prob[:i]) for i in range(0, len(prob))]) / sum(prob)
        df = pd.DataFrame({"t": x, "cumulative": cumulative})
        df["$\\alpha$"] = str(alpha)
        df_list.append(df)

    df = pd.concat(df_list).reset_index(drop=True)

    tau = 1000
    x = np.linspace(1, 20, num=1000)

    df_list = []
    for alpha in [0.5, 0.8, 1.5]:
        prob = alpha / (x ** (1 + alpha))
        cumulative = np.array([sum(prob[:i]) for i in range(0, len(prob))]) / sum(prob)
        df = pd.DataFrame({"t": x * tau, "cumulative": cumulative})
        df["$\\alpha$"] = str(alpha)
        df_list.append(df)

    df = pd.concat(df_list).reset_index(drop=True)

    sns.lineplot(data=df, x="t", y="cumulative", hue="$\\alpha$")
    plt.show()


if __name__ == "__main__":
    folder = "../production/pareto/"
    main(folder)
