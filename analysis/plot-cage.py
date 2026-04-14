# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pandas>=3.0.1",
#     "seaborn>=0.13.2",
# ]
# ///

import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set(font_scale=2)

from utils import get_msd_dataframe, add_pound_key, plateau_from_cage_size


def main(folder: str) -> None:
    files = glob.glob(f"{folder}/**/*.csv")

    df = pd.concat([pd.read_csv(f) for f in files]).reset_index(drop=True)

    # assert len(df) == len(files) * 45

    id_columns = ["time"]

    df_msd = get_msd_dataframe(df, id_columns)

    cage_size = 0.2

    # Save in wide format
    df_msd.pivot(index="time", values="MSD", columns="Definition").reset_index().to_csv(
        "cage_msd.csv", sep=" ", index=False
    )
    add_pound_key("cage_msd.csv")

    ax = sns.lineplot(data=df_msd, x="time", y="MSD", hue="Definition", linewidth=4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(y=plateau_from_cage_size(cage_size), color="grey", linestyle="dashed")
    plt.title("Cage")
    plt.show()


if __name__ == "__main__":
    folder = "../production/cage/"
    main(folder)
