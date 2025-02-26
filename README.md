# yukon-CK-ResDoc
Code to analyse yukon chinook data and generate corresponding CSAS Research Document

## To reproduce the analysis:
 - Fork or clone this repository  
 - Run the [code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/SR_all_CUs.R) to fit the spawner-recruit models, making sure you [toggle](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/SR_all_CUs.R#L12) `refit <- TRUE` on. This will also generate diagnostic and plots of the fits.  
 - Run the [forward simulation code](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/fwd_sims.R)
 - etc. as repo develops.

## For those who don't want to reproduce the analysis,  
you can view current [figures](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/tree/main/analysis/plots), [benchmarks](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/data/generated/bench_par_table.csv), and [data](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/tree/main/analysis/data/raw). This is subject to change as the analysis is updated. Once the CSAS Research Document is published, a final version of thie repository will be published on Zenodo. 
