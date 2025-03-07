bench.par1 <- read.csv(here("analysis/data/generated/bench_par_table.csv"))
bench.par2 <- read.csv(here("analysis/data/generated/brendan_bench_par_table.csv"))
bench.par3 <- read.csv(here("analysis/data/generated/hannah_bench_par_table.csv"))

#Want all these to be TRUE
all(bench.par1$CU == bench.par3$CU)
all(bench.par1$bench.par == bench.par3$bench.par)
all(bench.par1$mean == bench.par3$mean)

cbind(bench.par1$CU, bench.par3$CU)

#look at difference in means
cbind(bench.par1$CU, bench.par1$bench.par,  (bench.par1$mean - bench.par2$mean))
