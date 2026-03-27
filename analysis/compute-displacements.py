# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "h5py>=3.15.1",
#     "seaborn>=0.13.2",
# ]
# ///

import h5py
import numpy as np
import dataclasses
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


@dataclasses.dataclass
class SimulationResult:
    n_walker: int
    omegas: np.ndarray
    timesteps: np.ndarray
    H: float
    rate: float
    dt: float


def read_omegas(path: str) -> SimulationResult:
    f = h5py.File(path, "r")

    params = f["Params"]
    timesteps = params["scheduler"][:]

    omegas = []
    for time in timesteps:
        omega = np.array(f[f"TimeSteps/ExactRotation/{time}"])
        omega = np.array([list(o[0]) for o in omega])
        omegas.append(omega)

    n_walker = np.array(params["walkers"])
    H = np.array(params["H"])
    rate = np.array(params["rate"])
    dt = np.array(params["dt"])

    f.close()

    return SimulationResult(
        int(n_walker), np.array(omegas), timesteps, float(H), float(rate), float(dt)
    )


def type_from_h_rate(H: float, rate: float) -> str:
    if H == 0:
        if rate != 0:
            raise ValueError("Escape without cage??")
        return "Brownian"

    return "Cage" if rate == 0 else "Escape"


def main() -> None:
    files_to_process = glob.glob("examples/**/*.jld2", recursive=True)

    df_list = []
    for path in files_to_process:
        results = read_omegas(path)

        # We just look at the mean square norm of omega
        # omega is (time, walker, xyz), so we take the square norm on xyz, average over walkers
        omega_sq_norm = np.mean(np.linalg.norm(results.omegas, axis=2) ** 2, axis=1)
        df = pd.DataFrame(
            {"t": t * results.dt, "omega_sq_norm": o}
            for t, o in zip(results.timesteps, omega_sq_norm)
        )
        df["H"] = results.H
        df["rate"] = results.rate
        df["type"] = type_from_h_rate(results.H, results.rate)
        df_list.append(df)

    df = pd.concat(df_list).reset_index(drop=True)

    ax = sns.lineplot(data=df.query("type == 'Brownian'"), x="t", y="omega_sq_norm")
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.show()


if __name__ == "__main__":
    main()
