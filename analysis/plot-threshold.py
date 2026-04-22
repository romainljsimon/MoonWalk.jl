# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "seaborn>=0.13.2",
# ]
# ///

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set(font_scale=2)


def main(filepath: str) -> None:
    df = pd.read_csv(filepath)

    counter = 0
    jumps = [0]
    for i in range(1, len(df)):
        if abs(df["theta_local"][i] - df["theta_local"][i - 1]) > 1:
            counter += 1
        jumps.append(counter)
    df["nb_of_jump"] = [str(x) for x in jumps]

    threshold = 2

    def set_lines(ax, set_lims):
        # ax.set_ylabel("$\\theta$")
        ax.set_ylabel("")
        ax.axhline(y=threshold, linestyle="dashed", color="grey")
        ax.axhline(y=np.pi, linestyle="dashed", color="grey")
        if set_lims:
            ax.set_ylim([0, np.pi * 1.01])

    fig, ax = plt.subplots(3)

    sns.lineplot(data=df, x="time", y="theta", ax=ax[0])
    set_lines(ax[0], True)

    sns.lineplot(
        data=df,
        x="time",
        y="theta_local",
        hue="nb_of_jump",
        palette=sns.color_palette("tab10"),
        legend=False,
        ax=ax[1],
    )
    set_lines(ax[1], True)

    sns.lineplot(data=df, x="time", y="theta_rebuilt", ax=ax[2])
    set_lines(ax[2], False)

    for ax in fig.get_axes():
        ax.label_outer()

    plt.show()


if __name__ == "__main__":
    filepath = "../production/angles.csv"
    main(filepath)
