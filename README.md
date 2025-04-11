# yukon-CK-ResDoc
Code to analyze yukon chinook data and generate corresponding CSAS Research Document. Currently a work in progress.

## To reproduce the analysis:
 - Fork or clone this repository
 - Make sure all software [dependencies](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/790691908819bbc86733b1c1e72aea8dec73c0a5/analysis/R/load.R) are installed  
 - Run the [run-reconstruction model fitting code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/run-reconstructions/RR.R) to fit the run-reconstructions models, and generate time series of spawners, harvest and age, sex and size composition.
 - Run the [spawner-recruit model fitting code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/SR_fit.R) to fit the spawner-recruit models, making sure you [toggle](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/SR_fit.R#L8) `refit <- TRUE` on.
 - Run the [egg mass spawner-recruit model fitting code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/732b656ca908b2be3552e173b73f537cb82ef376/analysis/R/SR_fit_eggs.R) and [generate demographic reference points](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/732b656ca908b2be3552e173b73f537cb82ef376/analysis/R/demographic_ref_points.R)
 - Run the [forward simulation code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/fwd_sim.R)
 - Run the [figure and benchmark generation code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/inference_figs.R) and additional [tradeoff analyses](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/732b656ca908b2be3552e173b73f537cb82ef376/analysis/R/tradeoffs.R)

## For those who don't want to reproduce the analysis,  
you can view current [figures](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/tree/main/analysis/plots), [benchmarks](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/data/generated/bench_par_table.csv), and [data](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/tree/main/analysis/data/raw). This is subject to change as the analysis is updated. Once the CSAS Research Document is published, a final version of the repository will be published on Zenodo.

The model validation supplement is available [here](https://pacific-salmon-assess.github.io/yukon-CK-ResDoc/analysis/R/SR-diagnostics.html). 
