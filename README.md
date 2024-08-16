# 2024-frequency-memory

This repository contains the code that was used to generate the figures in [this preprint](https://arxiv.org/abs/2405.05180). 


## Structure of the repository

The repository contains the following:

- The `make_data.ipynb` notebook allows to create the great majority of the data required for the figures. This data is outputted in the `data` folder.
- The `figs.ipynb` notebook allows to create the minimal plots (axis, graduations, lines, heatmap, points, etc.) from the data in the `data` folder. It outputs these in the  `figures` folder.
- The `all_figs.svg` contains all the assembled plots using the assets in `figures`. The final figures were exported from [Inkscape](https://inkscape.org/) in the `final_figs` folder.
- The `individual_scripts` folder contains the scripts to generate the figures taking more time to generate (such as higher resolution Arnold tongues, as in Figure 6a). The data generated in this folder is used by the `make_data.ipynb` notebook.
- The `Manifest.toml` and `Project.toml` keep track of the versions of the packages used in to make the figures and the simulations.