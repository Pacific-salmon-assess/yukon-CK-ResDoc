# fit SR models for all CUs, run diagnostics, make initial plots
library(here)
source(here("analysis/R/data_functions.R"))

bench.posts <- readRDS(here("analysis/data/generated/benchmark_posteriors.rds"))
bench.long <- pivot_longer(bench.posts, cols = c(Sgen, Smsy, S.recent), names_to = "par") |>
  select(-Umsy, - Seq) |>
  arrange(CU, par, value) |>
  filter(value <= 10000) #hack to cut off fat tails to help with density visualization, also an IUCN cutoff... 

for(i in unique(sp_har$cu)){
  sub <- filter(bench.posts, CU == i)
  
  p <- ggplot(sub) +
    geom_density(aes(Smsy), color = "forestgreen", fill = "forestgreen", alpha = 0.2) +
    geom_density(aes(Sgen), color = "darkred", fill = "darkred", alpha = 0.2) +
    geom_density(aes(S.recent), fill = "grey", alpha = 0.5) +
    geom_vline(xintercept = 1500) +
    coord_cartesian(xlim = c(0, quantile(sub$Smsy, 0.99))) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    labs(x = "Spawners", y = "Posterior density", 
         title = paste(i, "spawner density and benchmarks"))
  print(p)
  my.ggsave(here("analysis/plots/CU_detail", paste0("status_", i, ".PNG")))
  #ggplot(sub) +
    #geom_density(aes(S.recent), fill = "grey", alpha = 0.2)
} 

ggplot(bench.long, aes(value/1000, fill = par, color = par)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 1.5) +
  facet_wrap(~CU, scales = "free_y") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("black", "darkred", "forestgreen"), 
                    aesthetics = c("fill", "color"), 
                    labels = c(expression(italic(S[recent])), expression(italic(S[gen])), 
                               expression(italic(paste("80% ",S)[MSY])))) +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title=element_blank()) +
  labs(x = "Spawners (thousands)", y = "Posterior density", 
       title = "Recent spawners relative to benchmarks and 1500 cutoff")
my.ggsave(here("analysis/plots/status.PNG"))
