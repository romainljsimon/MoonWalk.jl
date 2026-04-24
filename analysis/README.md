# Analysis scripts

    This folder contains python scripts used to postprocess the trajectories. All the production scripts (in the production subfolder) must be run before the analysis is carried out.


These scripts are meant to be used with [uv](https://docs.astral.sh/uv/). For example:
```bash
uv run -s plot-brownian.py
```
Alternatively, you can use a different python package manager, and install the dependencies manually. Each script has its own dependencies defined at the top of the file (seaborn, scipy, and plotly should be enough to run all the scripts).


Running these scripts will create some plots, which are similar to the figures in the publication (but not identical). They will also create some csv files with the post-processed data (typically MSD curves), which are exactly the data used in the publication's figures.
