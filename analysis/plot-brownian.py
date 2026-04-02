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

from utils import METHODS, get_rmsd_dataframe


def main(folder: str) -> None:
    files = glob.glob(f"{folder}/**/*.csv")

    df = pd.concat([pd.read_csv(f) for f in files]).reset_index(drop=True)

    assert len(df) == len(files) * 100

    id_columns = ["time"]

    df_rmsd = get_rmsd_dataframe(df, id_columns)

    # Theoretical value
    # Jump at times between [0, T], of an amplitude between [0, A]
    # RMSD is sqrt(A**2 * t / T**2)
    T = 1
    A = 0.1 / 2
    x = np.logspace(1, 4, 1000)
    D = 2 * A**2 / T
    y = np.sqrt(D * x)

    ax = sns.lineplot(data=df_rmsd, x="time", y="RMSD", hue="Definition", linewidth=4)
    ax.plot(x, y, linestyle="dashed", color="grey")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(y=np.sqrt((np.pi**2 + 6) / 3), color="grey", linestyle="dashed")
    plt.title("Brownian diffusion")
    plt.show()


if __name__ == "__main__":
    folder = "../production/brownian/"
    main(folder)
