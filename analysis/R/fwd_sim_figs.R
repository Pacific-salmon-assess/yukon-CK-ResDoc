
# forward simulations ----

## reference vs robustness productivity ----

AR1.par.posts <- read.csv(here("analysis/data/generated/AR1_posteriors.csv"))
TVA.par.posts <- read.csv(here("analysis/data/generated/TVA_posteriors.csv"))

TV.pp.ref.long <- pivot_longer(TVA.par.posts, cols = c(30:34), names_to = "par") |>
  select(CU, par, value) |>
  mutate(scenario = "reference (most recent generation)")

TV.pp.rob2.long <- pivot_longer(TVA.par.posts, cols = c(1:34), names_to = "par") |>
  select(CU, par, value) |>
  mutate(scenario = "robustness (long-term average)")

TV.pp.rob.long <- pivot_longer(TVA.par.posts, cols = c(34), names_to = "par") |>
  select(CU, par, value) |>
  mutate(scenario = "robustness (most recent year)")

AR.pp.rob <- AR1.par.posts |>
  mutate(scenario = "stationary",
         par= 1,
         value = ln_a) |>
  select(CU, par, value, scenario)

alpha.posts <- rbind(TV.pp.ref.long, TV.pp.rob.long, AR.pp.rob,TV.pp.rob2.long)

# Alpha comparison
alpha.stats <- alpha.posts |>
  group_by(CU, scenario) |>
  summarise(mean.prod = mean(value)) |>
  filter(scenario != "stationary")

# Plot comparison of scenarios
alpha.posts |> filter(scenario != "stationary") |>
  mutate(CU_f = factor(CU, levels=CU_order)) |>
  mutate(scenario = str_to_sentence(scenario)) |>
  ggplot(aes(value, fill = scenario, color = scenario)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~CU_f, scales = "free_y", labeller=CU_labeller) +
  theme_sleek() +
  scale_x_continuous(limits = c(NA, 4)) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, lty=2, col="grey") +
  scale_colour_grey(aesthetics = c("colour", "fill"),start = 0.3, end = 0.6) +
  labs(y = "", x = expression(Log(alpha)), fill="Productivity Scenario", color="Productivity Scenario")

my.ggsave(here("analysis/plots/fwd-sim/OM-productivity-scenarios.PNG"))
ggsave(here("csasdown/figure/OM-productivity-scenarios.PNG"), height=600*2, width=777*2,
       units = "px", dpi=240)

# Plot comparison of TVA long-term average vs stationary alpha
alpha.posts |> filter(scenario %in% c("stationary", "robustness (long-term average)")) |>
  mutate(scenario2 = case_when(scenario == "robustness (long-term average)" ~ "TV model average",
                               scenario == "stationary" ~ "Stationary model",
                               .default = scenario)) |>
  mutate(CU_f = factor(CU, levels=CU_order)) |>
  ggplot(aes(value, fill = scenario2, color = scenario2)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~CU_f, scales = "free_y", labeller=CU_labeller) +
  theme_sleek() +
  scale_x_continuous(limits = c(NA, 4)) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, lty=2, col="grey") +
  scale_colour_grey(aesthetics = c("colour", "fill"),start = 0.3, end = 0.6) +
  labs(y = "", x = expression(Log(alpha)), fill="Productivity Scenario", color="Productivity Scenario")
my.ggsave(here("analysis/plots/fwd-sim/OM-long-term-scenarios.PNG"))



# Load both S-R and egg mass -based biol. benchmarks by CU
bench.par.table <- read.csv(here("analysis/data/generated/bench_par_table.csv"))
demo.bench <- read.csv(here("analysis/data/generated/demographic_parameters.csv"))
# Make a dataframe of benchmarks used in sims
all.bench <- bench.par.table |>
  filter(bench.par == "Smsr") |>
  mutate(lwr_bench = 0.2*X50.,
         upr_bench = 0.4*X50.) |>
  left_join(filter(demo.bench, par=="Smsr", period=="recent"), by="CU") |>
  rename(rt = median) |>
  select(c("CU", "upr_bench", "lwr_bench", "rt"))


# Generate plots for which set of fwd simulations?
sceanarios <- list("TVA", "TVA2", "TVA3") # "AR1") # Can omit one to avoid re-generating figures

for(k in sceanarios) { # generate Fwd-sim figures for reference set (TVA) & robustness set (AR1)

  if(k == "TVA"){
    S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd_TVA.csv"))
    H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd_TVA.csv"))
    perf.metrics <- read.csv(here("analysis/data/generated/perf_metrics_TVA.csv")) |>
      pivot_longer(-c(1,12), names_to = "metric") |>
      pivot_wider(names_from=prob, values_from=value)
    Sig.R <- read.csv(here("analysis/data/generated/simulations/var_covar_TVA.csv"),
                      row.names = 1)
    spwn.obs = TV.spwn; harv.obs = TV.harv
  } else if (k == "TVA2"){
    S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd_TVA2.csv"))
    H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd_TVA2.csv"))
    perf.metrics <- read.csv(here("analysis/data/generated/perf_metrics_TVA2.csv")) |>
      pivot_longer(-c(1,12), names_to = "metric") |>
      pivot_wider(names_from=prob, values_from=value)
    Sig.R <- read.csv(here("analysis/data/generated/simulations/var_covar_TVA2.csv"),
                      row.names = 1)
    spwn.obs = TV.spwn; harv.obs = TV.harv
  } else if (k == "TVA3"){
    S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd_TVA3.csv"))
    H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd_TVA3.csv"))
    perf.metrics <- read.csv(here("analysis/data/generated/perf_metrics_TVA3.csv")) |>
      pivot_longer(-c(1,12), names_to = "metric") |>
      pivot_wider(names_from=prob, values_from=value)
    Sig.R <- read.csv(here("analysis/data/generated/simulations/var_covar_TVA3.csv"),
                      row.names = 1)
    spwn.obs = TV.spwn; harv.obs = TV.harv
  } else if (k == "AR1"){
    S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd_AR1.csv"))
    H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd_AR1.csv"))
    perf.metrics <- read.csv(here("analysis/data/generated/perf_metrics_AR1.csv")) |>
      pivot_longer(-c(1,12), names_to = "metric") |>
      pivot_wider(names_from=prob, values_from=value)
    Sig.R <- read.csv(here("analysis/data/generated/simulations/var_covar_AR1.csv"),
                      row.names = 1)
    spwn.obs = AR1.spwn; harv.obs = AR1.harv

  }

  # Variables as factors for plotting
  S.fwd$CU_f <- factor(S.fwd$CU, levels = CU_order)
  H.fwd$CU_f <- factor(H.fwd$CU, levels = CU_order)
  all.bench$CU_f <- factor(all.bench$CU, levels = CU_order)

  # assign HCRs to groups for subsetting harvest scenarios / HCRs
  ER_seq <- seq(5, 100, 5) # Must match ER_seq in "fwd_sim.R"
  HCRs <- c("no.fishing", "moratorium", "IMEG", "moratorium.cap", "IMEG.cap", "PA.alternative", paste0("fixed.ER.", ER_seq))
  HCR_grps <- list(base = c("no.fishing", "fixed.ER.60", "IMEG"),
                   moratorium = c("no.fishing", "moratorium", "moratorium.cap", "PA.alternative"),
                   IMEG = c("no.fishing",
                            unique(S.fwd$HCR[grepl("*IMEG*", S.fwd$HCR)])),
                   simple = c("no.fishing", "IMEG", "moratorium"),
                   fixed = unique(S.fwd$HCR[grepl("fixed.ER", S.fwd$HCR)]))
  # colours
  HCR_cols <- c("#B07300", "purple3", "grey25", "#CCA000", "#FEE106",  "#0F8A2E", "#3638A5")
  names(HCR_cols) <- c("PA Alternative", "Fixed ER 60", "No fishing", "Moratorium", "Moratorium cap", "IMEG", "IMEG cap")
  HCR_lookup <- data.frame(HCR=c(HCRs[1:6], "fixed.ER.60"), HCR_name = names(HCR_cols)[c(3:4,6,5,7,1,2)])
  HCR_lookup

  ## Spawners projection ----

  for(i in 1:length(HCR_grps[1:4])) { # don't make this fig for all fixed exp rates

    S.fwd |> filter(HCR %in% HCR_grps[[i]]) |>
      left_join(HCR_lookup, by="HCR") |>
      ggplot() +
      # Observations:
      geom_ribbon(data = filter(spwn.obs, year >= max(spwn.obs$year)-7),
                  aes(ymin = S.25/1000, ymax = S.75/1000,
                      x= year), #offset to return year
                  fill = "grey", color = "grey") +
      geom_line(data = filter(spwn.obs, year >= max(spwn.obs$year)-7),
                aes(y=S.50/1000, x= year), color = "black") +
      # Projections:
      geom_ribbon(aes(ymin = S.25/1000, ymax = S.75/1000, x = year, color=HCR_name, fill = HCR_name),
                  alpha = 0.2) +
      geom_line(aes(year, S.50/1000, color = HCR_name), lwd=1) +
      # Benchmarks:
      geom_hline(data=all.bench,
                 aes(yintercept = rt/1000), lty=2,
                 color = "pink3") +
      geom_hline(data=all.bench,
                 aes(yintercept = upr_bench/1000), lty=2,
                 color = "forestgreen") +
      geom_hline(data=all.bench,
                 aes(yintercept = lwr_bench/1000), lty=2,
                 color = "darkred") +
      scale_linetype_manual(values=2, guide = "legend") +
      facet_wrap(~CU_f, scales = "free_y", labeller=CU_labeller) +
      scale_x_continuous(expand = expansion(mult = c(0, .01))) +
      labs(y = "Spawners (000s)", x="Year", col="", fill="") +
      theme_sleek() +
      theme(legend.position = "bottom",
            strip.text = element_text(size=10),
            legend.text = element_text(size=10)) +
      scale_color_manual(values=HCR_cols, aesthetics = c("fill", "color"))

    my.ggsave(here(paste("analysis/plots/fwd-sim/S-fwd", names(HCR_grps[i]), "grp", paste0(k, ".PNG"), sep="_")))
    if(names(HCR_grps[i])=="simple" && k=="TVA"){
      ggsave(here(paste("csasdown/figure/S-fwd", names(HCR_grps[i]), "grp", paste0(k, ".PNG"), sep="_")),
             height=650*2, width=810*2, units="px", dpi=240)
    }
  }


  ## Harvest projection ----
  for(i in 1:length(HCR_grps[1:4])) { # don't make this fig for all fixed exp rates
    H.fwd |> filter(HCR %in% HCR_grps[[i]], HCR != "no.fishing") |>
      left_join(HCR_lookup, by="HCR") |>
      ggplot() +
      # Observations:
      geom_ribbon(data = filter(harv.obs, year >= max(harv.obs$year)-7),
                  aes(ymin = H.25/1000, ymax = H.75/1000,
                      x= year), #offset to return year
                  fill = "grey", color = "grey") +
      geom_line(data = filter(harv.obs, year >= max(harv.obs$year)-7),
                aes(y=H.50/1000, x= year), color = "black") +
      # Projections:
      geom_ribbon(aes(ymin = H.25/1000, ymax = H.75/1000, x = year, color=HCR_name, fill = HCR_name),
                  alpha = 0.2) +
      geom_line(aes(year, H.50/1000, color = HCR_name), lwd=1) +
      facet_wrap(~CU_f, scales = "free_y", labeller=CU_labeller) +
      scale_x_continuous(expand = expansion(mult = c(0, .01))) +
      labs(y = "Harvest (000s)", col="", fill="") +
      theme_sleek() +
      theme(legend.position = "bottom") +
      scale_color_manual(values=HCR_cols, aesthetics = c("fill", "color"))

    my.ggsave(here(paste("analysis/plots/fwd-sim/H-fwd", names(HCR_grps[i]), "grp", paste0(k, ".PNG"), sep="_")))
  }


  ## -- performance metrics and status multipanel (all HCR excl. fixed ER, all PMs) ----
  #metrics
  perf.metrics <- perf.metrics |> left_join(HCR_lookup, by="HCR") |>
    mutate(HCR_name = factor(HCR_name, levels=HCR_lookup$HCR_name))

  pm_plot <- perf.metrics |>
    filter(!(HCR %in% HCR_grps[["fixed"]])) |>
    filter(!(HCR %in% "realistic")) |>
    filter(!(metric %in% c("n.above.upr", "n.between.bench", "n.below.lwr", "n.above.reb", "n.extinct"))) |>
    mutate(metric_name = case_when(metric == "ER" ~ "Exploitation rate",
                                   metric == "cdn.harvest" ~ "Canadian harvest",
                                   metric == "pr.closed" ~ "Prop. years fishery closed",
                                   .default = str_to_sentence(metric)))|>
    mutate(metric_name = factor(metric_name, levels=c("Escapement", "Harvest", "Canadian harvest", "Prop. years fishery closed","Exploitation rate"))) |>
    ggplot() +
    geom_col(aes(x=HCR_name, y = median, fill=HCR_name)) +
    geom_segment(aes(x=HCR_name,
                     xend=HCR_name,
                     y=q_25, yend=q_75), col="grey30") +
    scale_fill_manual(values=HCR_cols) +
    facet_wrap(~metric_name, scales = "free_y", nrow=4) +
    theme_sleek() +
    scale_y_continuous() +
    theme(legend.position = c(0.8,0.15),
          axis.text.x = element_blank(),
          legend.title = element_blank()) +
    guides(fill=guide_legend(ncol=2, byrow=T)) +
    labs(x="", y="")

  # status
  perf.status <- perf.metrics |>
    filter(metric %in% c("n.above.upr", "n.between.bench", "n.below.lwr", "n.above.reb", "n.extinct")) |>
    mutate(status = factor(gsub("^n ", "", gsub("\\.", " ", metric)),
                           levels=c("above reb", "above upr", "between bench", "below lwr", "extinct")))

  status_plot <- perf.status |>
    filter(!(HCR %in% HCR_grps[["fixed"]])) |>
    ggplot(aes(x = HCR_name, y = mean, fill = status)) +
    geom_col() +
    scale_fill_discrete(type = c("pink3", "forestgreen", "darkorange", "darkred", "black")) +
    scale_y_continuous(breaks = c(2,4,6,8)) +
    labs(x="Harvest control rule", y = "Average # of CUs \n (1000 simulations)", fill="Status") +
    theme_sleek()

  cowplot::plot_grid(pm_plot, status_plot, nrow=2, labels="auto", rel_heights = c(1.5,1))

  my.ggsave(here(paste0("analysis/plots/fwd-sim/perf_metrics_status_", k, ".PNG")))
  ggsave(here(paste0("csasdown/figure/perf_metrics_status_", k, ".PNG")), height=900*2,
         width=800*2, units="px", dpi=240,bg = "white")


  ## Identify LRP-F --
  perf.metrics %>%
    filter(HCR %in% HCR_grps[["fixed"]]) %>%
    group_by(HCR) %>%
    mutate(ER = as.numeric(gsub("\\D", "", HCR))) %>%
    filter(metric == "n.below.LSR") %>%
    #mutate(prec = if_else(median <= 1, 1, 0)) %>%
    ggplot(aes(x=ER, y=mean)) +
    geom_point()


  ## fixed ER trade-off multipanel ----

  spwn_v_ER <- S.fwd %>% filter(HCR %in% HCR_grps[["fixed"]]) %>%
    group_by(HCR, CU_f) %>%
    summarize(mean_spwn = mean(S.50)) %>%
    mutate(ER = as.numeric(gsub("\\D", "", HCR))) %>%
    ggplot() +
    geom_point(aes(y=ER, x=mean_spwn/1000, col=CU_f), shape='circle', size=2, alpha=0.7) +
    scale_colour_viridis_d() + scale_y_continuous(breaks = seq(0,100,20)) +
    theme_sleek() +
    theme(legend.position="none",
          axis.title = element_text(size=12)) +
    labs(x="Spawners (000s)", y="Exploitation rate", col="Conservation Unit")

  harv_v_ER <- H.fwd |> filter(HCR %in% HCR_grps[["fixed"]]) |>
    group_by(HCR, CU_f) |>
    summarize(mean_harv = mean(H.50)) |>
    mutate(ER = as.numeric(gsub("\\D", "", HCR))) |>
    left_join(CU_name_lookup) |>
    ggplot() +
    geom_point(aes(y=ER, x=mean_harv/1000, col=CU_pretty), shape='circle', size=2, alpha=0.7) +
    scale_colour_viridis_d() + scale_y_continuous(breaks = seq(0,100,20)) +
    theme_sleek() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.text = element_text(size=10),
          axis.title = element_text(size=12)) +
    labs(x="Harvest (000s)", y="Exploitation rate", col="Conservation Unit")

  status_ER <- perf.status %>% filter(HCR %in% HCR_grps[["fixed"]]) %>%
    mutate(ER = as.numeric(gsub("\\D", "", HCR))) %>%
    filter(ER != 100) %>%
    ggplot(aes(x=mean, y=factor(ER), fill=status)) +
    geom_col() +
    geom_vline(data=data.frame(x=seq(1:9)), aes(xintercept=x), col="white", linewidth=0.05) +
    scale_fill_discrete(type = c("pink3","forestgreen", "darkorange", "darkred", "black")) +
    scale_x_continuous(breaks = c(2,4,6,8)) +
    scale_y_discrete(breaks = seq(0,90,20)) +
    labs(y = "Exploitation Rate", fill="", x="Average # of CUs (over 1000 simulations)") +
    theme_sleek() +
    theme(legend.margin = margin(l=20, r=20),
          legend.text = element_text(size=10),
          axis.title = element_text(size=12))

  b <- cowplot::plot_grid(spwn_v_ER, harv_v_ER, nrow=1, rel_widths=c(0.55,1), labels=c("b","c"), label_x = c(0,-0.015))
  cowplot::plot_grid(status_ER, b, nrow=2, rel_heights=c(1,1), labels=c(NULL,"a"))

  my.ggsave(here(paste0("analysis/plots/fwd-sim/fixed_ER_tradeoffs_", k, ".PNG")))
  ggsave(here(paste0("csasdown/figure/fixed_ER_tradeoffs_", k, ".PNG")), height=600*2,
         width=800*2, units="px", dpi=240)


  ## visualize HCRs ----
  HCR_order <- c("Moratorium", "IMEG", "Moratorium cap", "IMEG cap", "PA Alternative")
  out <- visualize_HCR(HCRs=HCRs[2:6]) # get simulated HRs
  out <- left_join(out, HCR_lookup, by="HCR")
  out$HCR_name <- factor(out$HCR_name, levels=HCR_order)

  ggplot(out) + geom_line(aes(x=run_size/1000, y=HR*100, col=HCR_name), linewidth=0.75) +
    scale_colour_manual(values=HCR_cols, guide="none") +
    scale_fill_manual(values="grey70", guide="legend") +
    facet_wrap(~HCR_name) +
    labs(x="Run Size (000s)", y="Harvest Rate (%)") +
    theme_minimal() + theme(strip.text = element_text(size=10),
                            axis.title = element_text(size=10)) +
    lims(x=c(0,200)) +
    scale_y_continuous(breaks=seq(0,100,20), limits=c(0,100))

  my.ggsave(here("analysis/plots/fwd-sim/HCR_visualize.PNG"))
  my.ggsave(here("csasdown/figure/HCR_visualize.PNG"), height=4,
            width=6)


  ## covariance plots ----

  Sig.R.order <- Sig.R[c(4,8,6,2,5,3,1,7,9),c(4,8,6,2,5,3,1,7,9)]

  ggcorrplot(Sig.R.order, hc.order = TRUE, type = "lower",
             outline.col = "white",
             lab=TRUE) +
    scale_fill_gradient2(limit=c(0,max(Sig.R[row(Sig.R) != col(Sig.R)])),
                         high="red") # this scales the colours for each plot

  my.ggsave(here(paste0("analysis/plots/fwd-sim/recruit-corr-matrix_", k, ".PNG")))

} # end k loop
