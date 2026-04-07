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
    get_msd_dataframe,
    plateau_from_cage_size,
    get_diffusion_coefficient,
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

    id_columns = ["time", "alpha"]

    df_msd = get_msd_dataframe(df, id_columns)

    cage_size = 0.2
    plateau = plateau_from_cage_size(cage_size)

    df_msd["$\\alpha$"] = df_msd["alpha"].apply(str)

    # All alphas, unbounded
    # ax = sns.lineplot(
    #    data=df_msd.query("Definition == 'Unbounded'"),
    #    x="time",
    #    y="MSD",
    #    hue="$\\alpha$",
    #    linewidth=4,
    # )
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # ax.axhline(y=np.sqrt(plateau), color="black", linestyle="dashed")
    # ax.plot(x, y, linestyle="dashed", color="grey")
    # plt.title("Pareto")
    # plt.show()

    # one alpha, look at all definitions
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
