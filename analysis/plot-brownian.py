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

from utils import METHODS, get_msd_dataframe


def main(folder: str) -> None:
    files = glob.glob(f"{folder}/**/*.csv")

    df = pd.concat([pd.read_csv(f) for f in files]).reset_index(drop=True)

    assert len(df) == len(files) * 100

    id_columns = ["time"]

    df_msd = get_msd_dataframe(df, id_columns)

    # Save in wide format
    df_msd.pivot(index="time", values="MSD", columns="Definition").reset_index().to_csv(
        "brownian_msd.csv", sep=" ", index=False
    )

    # Theoretical value
    # Jump at times of mean T, of an amplitude between [0, A]
    # MSD is sqrt(A**2 * t / T)
    T = 1 / 1  # exponential of lambda: mean is 1/lambda
    A = 0.1 / 2
    x = np.logspace(-1, 4, 1000)
    D = A**2 / T
    y = D * x

    ax = sns.lineplot(data=df_msd, x="time", y="MSD", hue="Definition", linewidth=4)
    ax.plot(x, y, linestyle="dashed", color="grey")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(y=(np.pi**2 + 6) / 3, color="grey", linestyle="dashed")
    plt.title("Brownian diffusion")
    plt.show()


if __name__ == "__main__":
    folder = "../production/brownian/"
    main(folder)
