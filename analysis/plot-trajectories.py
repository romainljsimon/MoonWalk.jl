# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pandas>=3.0.1",
#     "plotly>=6.6.0",
# ]
# ///

import numpy as np
import pandas as pd
import plotly.graph_objects as go


def plot_trajectory(filepath):
    df = pd.read_csv(filepath)

    u = np.linspace(0, 2 * np.pi, 500)
    v = np.linspace(0, np.pi, 500)

    R = 0.99  # Radius of the sphere
    xs = R * np.outer(np.cos(u), np.sin(v))
    ys = R * np.outer(np.sin(u), np.sin(v))
    zs = R * np.outer(np.ones_like(u), np.cos(v))

    fig = go.Figure()
    # Sphere
    fig.add_trace(
        go.Surface(
            x=xs,
            y=ys,
            z=zs,
            # colorscale=[[0.0, "#f7fbff"], [0.5, "#c6dbef"], [1.0, "#9ecae1"]],
            colorscale=[[0.0, "#9ecae1"], [1.0, "#9ecae1"]],
            # colorscale="Blues",
            showscale=False,
            opacity=1.0,
            lighting=dict(
                ambient=0.6,
                diffuse=0.7,
                specular=0.05,  # Higher means stronger light dot
                roughness=0.4,  # Higher means larger dot
                fresnel=0.05,
            ),
            lightposition=dict(
                x=0,
                y=0,
                z=5,
            ),
        )
    )
    # Trajectory
    fig.add_trace(
        go.Scatter3d(
            x=df["x"],
            y=df["y"],
            z=df["z"],
            mode="lines",
            line=dict(
                color=df["time"],
                width=6,
                colorscale="Inferno",  # Viridis, Plasma, Inferno, Cividis
            ),
        )
    )

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
        )
    )

    fig.show()


def main() -> None:
    # filepath = "../production/trajectory/brownian/positions.csv"
    # plot_trajectory(filepath)
    # filepath = "../production/trajectory/cage/positions.csv"
    # plot_trajectory(filepath)
    filepath = "../production/trajectory/escape/positions.csv"
    plot_trajectory(filepath)


if __name__ == "__main__":
    main()
