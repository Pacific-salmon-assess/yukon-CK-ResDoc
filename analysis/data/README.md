# Data sources  

- 📁`raw` includes raw data and data generated in models needed to create this analysis.  
    - 📁`subfolder_1`   
    - 📁`subfolder_2`   

- 📁`generated` data generated from the analyses in this repo  
| File | Description | Generated in |
|------|-------------|--------------|
|`bench_par_table.csv` | Spawner-Recruit model parameters and biological benchmarks with their associated uncertainties and some diagnostics |  [inference-figs.R](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/analysis/R/inference-figs.R) |
| `perf_metrics.csv`| Table of performance metrics for each harvest control rule evaluated in the the forward simulation | [fwd.sim.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/fwd_sim.R) | 
| `TV_refpts.csv` | Time-varying reference points calculated using the most recent generation's (i.e., 5 year) average estimate of productivity ***to be deleted?***| [TV_ref_points.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/TV_ref_points.R)|
| 📁`simulations` | folder containing files `H_fwd.csv` and `S_fwd.csv` which describe projections of catch (i.e., harvest) and spawners through time under different HCRs | [fwd.sim.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/fwd_sim.R) |
