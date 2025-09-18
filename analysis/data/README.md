# Data sources  

## 📁`raw` includes raw data and data generated in models needed to create this analysis.  

## 📁`generated` data are files generated from the analyses in this repository.  

| File | Description | Generated in |  
| ---- | ----------- | ------------ |  
| `brood_table_long.csv` | Brood table (spawners and recruits) for CUs | [inference-figs.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/inference-figs.R) |
| `bench_par_table.csv` | Spawner-Recruit model parameters and biological benchmarks with their associated uncertainties and some diagnostics |  [inference-figs.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/inference-figs.R) | 
| `demographic_parameters.csv` | Biological benchmarks derived from demographic (egg mass) models | [demographic_ref_points.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/demographic_ref_points.R) |
| `benchmark_posteriors.rds` | Full posteriors od spawner-Recruit model parameters and biological benchmarks used to propogate uncertainty through forward simulations |  [inference-figs.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/inference-figs.R) |  
| `perf_metrics_x.csv` | Table of performance metrics for each harvest control rule evaluated in the the forward simulation. End of file name ('x') is the productivity scenario forward simulations were conditioned on. | [fwd.sim.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/fwd_sim.R) |   
| `TV_refpts.csv` | Time-varying reference points calculated using the most recent generation's (i.e., 5 year) average estimate of productivity ***to be deleted?*** | [TV_ref_points.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/TV_ref_points.R) |  
| 📁`simulations` | folder containing files `H_fwd_x.csv` and `S_fwd_x.csv` which describe projections of catch (i.e., harvest) and spawners through time under different HCRs. End of file name ('x') is the productivity scenario forward simulations were conditioned on. | [fwd.sim.R](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/fwd_sim.R) |  


    Productivity scenarios:
        - TVA: Recent generational average of _alpha_ estimates from time-varying alpha (TVA) models ('reference scenario')
        - TVA2: Most recent year _alpha_ estimate from time-varying alpha (TVA) models ('worst case scenario')
        - TVA3: Long-term average of _alpha_ estimates from time-varying alpha (TVA) models ('best case scenario')
        - AR1: _alpha_ estimate from stationary, autoregressive (AR1) models (not used)
