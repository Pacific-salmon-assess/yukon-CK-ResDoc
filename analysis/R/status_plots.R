# fit SR models for all CUs, run diagnostics, make initial plots
library(here)
source(here("analysis/R/data_functions.R"))

bench.posts <- readRDS(here("analysis/data/generated/benchmark_posteriors.rds"))
#bench.long <- pivot_longer()

for(i in unique(sp_har$cu)){
  sub <- filter(bench.posts, CU == i)
  
  p <- ggplot(sub) +
    geom_density(aes(Smsy), color = "forestgreen", fill = "forestgreen", alpha = 0.2) +
    geom_density(aes(Sgen), color = "darkred", fill = "darkred", alpha = 0.2) +
    geom_density(aes(S.recent), fill = "grey", alpha = 0.2) +
    geom_vline(xintercept = 1500) +
    coord_cartesian(xlim = c(0, quantile(sub$Smsy, 0.99))) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    labs(x = "Spawners", y = "Posterior density", 
         title = paste(i, "spawner density and benchmarks"))
  print(p)
  my.ggsave(here("analysis/plots", paste0("status_", i, ".PNG")))
  #ggplot(sub) +
    #geom_density(aes(S.recent), fill = "grey", alpha = 0.2)
} 
