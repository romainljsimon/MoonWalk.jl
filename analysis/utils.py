import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

METHODS = ["ExactRotation", "Integral", "Unbounded"]


def get_rmsd_dataframe(df: pd.DataFrame, id_columns=list[str]) -> pd.DataFrame:
    for method in METHODS:
        df[f"{method}"] = (
            df[f"{method}_x"] ** 2 + df[f"{method}_y"] ** 2 + df[f"{method}_z"] ** 2
        )

    return (
        df[id_columns + METHODS]
        .melt(id_vars=id_columns, value_vars=METHODS, var_name="Definition")
        .groupby(id_columns + ["Definition"])
        .apply(lambda x: np.sqrt(np.mean(x)))
        .rename("RMSD")
        .reset_index()
    )


def plateau_from_cage_size(cage_size: float) -> float:
    return (
        -(3 * cage_size**2 - 6) * np.sin(cage_size)
        - 6 * cage_size * np.cos(cage_size)
        + cage_size**3
    ) / (3 * (cage_size - np.sin(cage_size)))


def get_diffusion_coefficient(
    ddf: pd.DataFrame, expected_value: float | None = None
) -> float:
    ddf["MSD"] = ddf["RMSD"] ** 2
    ddf["MSD_over_time"] = ddf["MSD"] / ddf["time"]
    D = np.mean(ddf["MSD_over_time"][-10:])

    # ax = sns.lineplot(data=ddf, x="time", y="MSD_over_time")
    # ax.set_xscale("log")
    # if expected_value:
    #    ax.axhline(y=expected_value, color="grey", linestyle="dashed")
    # ax.axhline(y=D, color="blue", linestyle="dashed")
    # plt.show()

    return D
