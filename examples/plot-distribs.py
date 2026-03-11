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
import scipy
import seaborn as sns

sns.set(font_scale=2)

methods = ["ExactRotation", "Integral", "Unbounded"]


def get_brownian_percentiles(percs):
    files = glob.glob("brownian/*/*.csv")
    df = pd.concat([pd.read_csv(f) for f in files]).reset_index(drop=True)
    time = 10000.0
    axes = ["Unbounded_x", "Unbounded_z", "Unbounded_z"]
    ddf = df[df["time"] == time]
    ddf = ddf[axes].melt(value_vars=axes, value_name="angle").drop(columns="variable")
    ddf["normalized_angle"] = ddf["angle"] / np.std(ddf["angle"])

    df_loc = pd.DataFrame(
        {
            "percentile": percs,
            "value": np.percentile(ddf["normalized_angle"], percs * 100),
        }
    )
    df_loc["Simulation"] = "Brownian"
    return df_loc


def get_pereto_percentiles(percs):
    files = glob.glob("pareto/*/*/*.csv")

    df_list = []
    for f in files:
        ddf = pd.read_csv(f)
        ddf["alpha"] = f.split("/")[1]
        df_list.append(ddf)

    df = pd.concat(df_list).reset_index(drop=True)

    time = 100000.0
    axes = ["Unbounded_x", "Unbounded_z", "Unbounded_z"]
    ddf = df[df["time"] == time]
    ddf = (
        ddf[["alpha"] + axes]
        .melt(id_vars="alpha", value_vars=axes, value_name="angle")
        .drop(columns="variable")
    )
    # Normalize by std
    ddf = ddf.merge(ddf.groupby("alpha").apply(np.std).rename("std").reset_index())
    ddf["normalized_angle"] = ddf["angle"] / ddf["std"]

    def get_percs(ddf, **kwargs):
        return pd.DataFrame(
            {
                "percentile": percs,
                "value": np.percentile(ddf["normalized_angle"], percs * 100),
            }
        )

    df_loc = ddf.groupby("alpha").apply(get_percs).reset_index().drop(columns="level_1")
    df_loc["Simulation"] = df_loc["alpha"].apply(lambda x: f"CTRW, alpha = {x}")
    return df_loc.drop(columns="alpha")


def main(folder: str) -> None:
    percs = np.linspace(1e-3, 1 - 1e-3, 100)

    df_brownian = get_brownian_percentiles(percs)
    df_pareto = get_pereto_percentiles(percs)

    df_perc = pd.concat([df_brownian, df_pareto]).reset_index(drop=True)

    df_perc = df_perc.merge(
        pd.DataFrame({"percentile": percs, "theoretical": scipy.stats.norm.ppf(percs)})
    )

    sns.lineplot(
        data=df_perc, x="theoretical", y="value", hue="Simulation", linewidth=4
    )
    plt.axline((0, 0), slope=1, linestyle="dashed", color="grey")
    plt.title("QQ-plot")
    plt.xlabel("Normal quantiles")
    plt.ylabel("Measured quantiles")
    plt.show()


if __name__ == "__main__":
    main()
