# yukon-CK-ResDoc
Code to analyze yukon chinook data and generate corresponding CSAS Research Document. This is currently a WIP.

## To reproduce the analysis:
 - Fork or clone this repository
 - Make sure all software **dependencies** (add) are installed  
 - Run the [spawner-recruit model fitting code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/SR_fit.R) to fit the spawner-recruit models, making sure you [toggle](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/SR_fit.R#L8) `refit <- TRUE` on. 
 - Run the [forward simulation code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/fwd_sims.R)
 - Run the [figure and benchmark generation code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/inference_figs.R)

## For those who don't want to reproduce the analysis,  
you can view current [figures](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/tree/main/analysis/plots), [benchmarks](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/data/generated/bench_par_table.csv), and [data](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/tree/main/analysis/data/raw). This is subject to change as the analysis is updated. Once the CSAS Research Document is published, a final version of thie repository will be published on Zenodo.

The model validation supplement is available [here](https://pacific-salmon-assess.github.io/yukon-CK-ResDoc/analysis/R/SR-diagnostics.html). 
