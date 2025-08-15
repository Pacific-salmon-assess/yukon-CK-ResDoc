## Exploratory code to determine FRP-L

spwn <- read.csv(here('analysis', 'data', 'generated', 'simulations',
                               paste0('S_fwd_', 'TVA3', '.csv')))

sr.bench <- read.csv(here("analysis/data/generated/bench_par_table.csv")) |>
  filter(bench.par == "Smsr") |>
  select(c("CU", "bench.par", "X50.")) |>
  mutate(lb_med = X50.*0.2)

# Calculate CUs below lower benchmark
spwn.sum <-
  spwn |> left_join(sr.bench[,c(1,4)], by=c("CU")) |>
  mutate(below = if_else(S.50 < lb_med, TRUE, FALSE)) |> # is esc<BM?
  group_by(HCR, year) |>
  summarize(num_below=sum(below), smu_esc=sum(S.50)) |> # how many CU<BM?
  filter(HCR %in% HCR_grps[["fixed"]])
View(spwn.sum)

# Plot # of CUs<BM vs SMU escapement
spwn.sum |>
  mutate(ER = as.numeric(gsub("\\D", "", HCR))) |>
  ggplot(aes(x=smu_esc, y=num_below, col=year)) + geom_point() +
  scale_colour_viridis() +
  scale_y_continuous(breaks=seq(0,9,1)) +
  lims(x=c(0,50000))

# Average SMU escapement when median of 1 CU in red zone (50% prob)
spwn.sum |> ungroup() |> filter(num_below == 1) |>
  summarize(min=min(smu_esc), max=max(smu_esc), mean=mean(smu_esc))

