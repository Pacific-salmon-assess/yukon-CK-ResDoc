bench.par1 <- read.csv(here("analysis/data/generated/bench_par_table.csv"))
bench.par2 <- read.csv(here("analysis/data/generated/brendan_bench_par_table.csv"))
bench.par3 <- read.csv(here("analysis/data/generated/hannah_bench_par_table.csv"))

#Want all these to be TRUE, run inside of all() to see WHICH dont match
all(bench.par1$CU == bench.par2$CU)
#if above is FALSE, which CUs are getting mixed up?
cbind(bench.par1$CU, bench.par2$CU)[which(bench.par1$CU != bench.par2$CU),]

all(bench.par1$bench.par == bench.par2$bench.par)
all(bench.par1$mean == bench.par2$mean) #betas match

cbind(bench.par1$CU, bench.par3$CU)

#look at difference in means
cbind(bench.par1$CU, bench.par1$bench.par,  (bench.par1$mean - bench.par2$mean))
cbind(bench.par1$CU, bench.par2$CU)

#is it just the order that is bessed up? sort means and compare
cbind(sort(bench.par1$mean), sort(bench.par2$mean))
