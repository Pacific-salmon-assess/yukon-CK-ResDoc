# yukon-CK-ResDoc
Code to reproduce analysis and documentation for: 
> Connors, B.M., O’Dell, A., Hunter, H., Glaser, D., Gill, J., Rossi, S., and Churchland, C. 2025. Stock status and biological and fishery consequences of alternative harvest and rebuilding actions for Yukon River Chinook salmon (*Oncorhynchus tshawytscha*). DFO Can. Sci. Advis. Sec. Res. Doc. 2025/nnn. iv + 130 p.

An in press version of the Research Document is available [here](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/dba18480b9a538ac8aec31df5903ccdb26448b88/Yukon-chinook-resdoc.4Nov2025.pdf). And presentation summarizing the work at regional peer-review meeting is available [here](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/Yukon%20Chinook%20CSAS.4-5June2025.for%20distribution.pdf).

## To reproduce the analysis:
 - Clone or fork this repository
 - Recreate the project environment ([`load_environment.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/load_environment.R)), then restart your R session.
 - Process the raw Genetic Stock Identification ([`gsi_wrangle.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/gsi_wrangle.R)) and age, sex, length ([`asl_wrangle.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/asl_wrangle.R)) datasets.
 - Run the run-reconstruction model fitting code ([`RR.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/run-reconstructions/RR.R)) to fit the multi-Conservation Unit border passage reconstructions model, and generate time series of spawners, harvest and age, sex, and size composition.
 - Run the spawner-recruit model fitting code ([`SR_fit.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/SR_fit.R)) to fit the spawner-recruit models, making sure you toggle `refit <- TRUE` on.
 - Run the egg mass spawner-recruit model fitting code ([`SR_fit_eggs.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/732b656ca908b2be3552e173b73f537cb82ef376/analysis/R/SR_fit_eggs.R)) and generate demographic reference points ([`demographic_ref_points.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/732b656ca908b2be3552e173b73f537cb82ef376/analysis/R/demographic_ref_points.R))
 - Run the forward simulation code ([`fwd_sim.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/fwd_sim.R))
 - Run the figure and benchmark generation code ([`inference_figs.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/inference_figs.R)) to generate key plots and tables for inference.
 - Run the forward simulation figure generation code ([`fwd_sim_figs.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/fwd_sim_figs.R))

## To reproduce the CSAS Research Document working paper:
 - Install [csasdown](https://github.com/pbs-assess/csasdown), and ensure you also have dependencies installed
 - Open [`index.Rmd`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/tree/main/csasdown) in RStudio and then click the knit button
 - The resulting pdf of the working paper will then be saved to a subfolder that is created called `book` 

## For those who don't want to reproduce the analysis:  
You can view current [figures](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/tree/main/analysis/plots), [benchmarks](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/data/generated/bench_par_table.csv), and [data](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/tree/main/analysis/data/raw). This is subject to change as the analysis is updated. Once the CSAS Research Document is published, a final version of the repository will be published on Zenodo.

The model validation supplement is available [here](https://pacific-salmon-assess.github.io/yukon-CK-ResDoc/analysis/R/SR-diagnostics.html). 
