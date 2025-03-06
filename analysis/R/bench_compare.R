bench.par1 <- read.csv(here("analysis/data/generated/bench_par_table.csv"))

bench.par2 <- read.csv(here("analysis/data/generated/hannah_bench_par_table.csv"))

cbind(bench.par1$CU, bench.par1$bench.par, (bench.par1$mean - bench.par2$mean))
