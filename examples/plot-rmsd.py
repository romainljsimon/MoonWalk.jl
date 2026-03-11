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

methods = ["ExactRotation", "Integral", "Unbounded"]


def main(folder: str) -> None:
    files = glob.glob(f"{folder}/**/*.csv")

    df = pd.concat([pd.read_csv(f) for f in files]).reset_index(drop=True)

    assert len(df) == len(files) * 45

    for method in methods:
        df[f"{method}"] = np.sqrt(
            df[f"{method}_x"] ** 2 + df[f"{method}_y"] ** 2 + df[f"{method}_z"] ** 2
        )

    df_rmsd = (
        df[["time"] + methods]
        .melt(id_vars="time", value_vars=methods, var_name="method")
        .groupby(["time", "method"])
        .apply(np.mean)
        .rename("RMSD")
        .reset_index()
    )

    ax = sns.lineplot(data=df_rmsd, x="time", y="RMSD", hue="method")
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} folder")
    main(sys.argv[1])
