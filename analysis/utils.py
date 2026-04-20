import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

METHODS = ["ExactRotation", "Integral", "Unbounded"]


def get_msd_dataframe(df: pd.DataFrame, id_columns=list[str]) -> pd.DataFrame:
    for method in METHODS:
        df[f"{method}"] = (
            df[f"{method}_x"] ** 2 + df[f"{method}_y"] ** 2 + df[f"{method}_z"] ** 2
        )

    return (
        df[id_columns + METHODS]
        .melt(id_vars=id_columns, value_vars=METHODS, var_name="Definition")
        .groupby(id_columns + ["Definition"])
        .apply(lambda x: np.mean(x))
        .rename("MSD")
        .reset_index()
    )


def plateau_from_cage_size(cage_size: float) -> float:
    return (
        -(3 * cage_size**2 - 6) * np.sin(cage_size)
        - 6 * cage_size * np.cos(cage_size)
        + cage_size**3
    ) / (3 * (cage_size - np.sin(cage_size)))


def get_diffusion_coefficient(
    ddf: pd.DataFrame,
    expected_value: float | None = None,
    plot: bool = False,
    title: str | None = None,
) -> float:
    ddf["MSD_over_time"] = ddf["MSD"] / ddf["time"]
    D = np.mean(ddf["MSD_over_time"][-7:]) / 3

    if plot:
        # Only get the last 3 decades
        df_to_plot = ddf[ddf["time"] >= ddf["time"].max() / 1000]
        ax = sns.lineplot(data=df_to_plot, x="time", y="MSD_over_time")
        ax.set_xscale("log")
        if expected_value:
            ax.axhline(y=expected_value, color="grey", linestyle="dashed")
        ax.axhline(y=D * 3, color="blue", linestyle="dashed")
        if title:
            plt.title(title)
        plt.ylabel("MSD / t")
        plt.show()

    return D


def add_pound_key(filepath: str) -> None:
    lines = open(filepath).readlines()
    lines[0] = "# " + lines[0]
    with open(filepath, "w") as f:
        [f.write(line) for line in lines]
