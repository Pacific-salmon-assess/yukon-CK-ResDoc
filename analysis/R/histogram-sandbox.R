bench.long <- pivot_longer(bench.posts, cols = c(Smsr.20, Smsr.40, S.recent), names_to = "par") |>
  arrange(CU, par, value) |>
  filter(value <= 10000) #hack to cut off fat tails to help with density visualization, also an IUCN cutoff... 

bs <- bench.long |>
  filter(CU == "Big.Salmon",
         value <= 6000)

n <- bench.long |>
  filter(CU == "Nordenskiold",
         value <= 2500)

s <- bench.long |>
  filter(CU == "Stewart",
         value <= 4000)

u <- bench.long |>
  filter(CU == "UpperYukonR.",
         value <= 4000)


w <- bench.long |>
  filter(CU == "Whiteandtribs.",
         value <= 5000)

t <- bench.long |>
  filter(CU == "YukonR.Teslinheadwaters",
         value <= 5000)

p <- bench.long |>
  filter(CU == "Pelly",
         value <= 9000)

m <- bench.long |>
  filter(CU == "MiddleYukonR.andtribs.",
         value <= 10000)

no <- bench.long |>
  filter(CU == "NorthernYukonR.andtribs.",
         value <= 10000)

custom.bench <- rbind(bs,n,s,u,w,t,p,m,no) |>
  mutate(CU_f = CU) |>
  left_join(CU_name_lookup, by="CU_f")

b <- ggplot(custom.bench |> filter(), aes(Smsr/1000, fill = CU_pretty, color = CU_pretty)) +
  geom_density(alpha = 0.3,bw=0.6) +
  theme(legend.position = "bottom") +
  labs(x = expression(italic(S[MSR])), y = "Posterior density") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_sleek()   +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title=element_blank(),
        legend.position="none") +
  scale_x_continuous(limits = c(0, 25))

c <- ggplot(custom.bench |> filter(), aes(Umsy, fill = CU_pretty, color = CU_pretty)) +
  geom_density(alpha = 0.3,bw=0.03) +
  theme(legend.position = "bottom") +
  labs(x = expression(italic(U[MSY])), y = "Posterior density") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_sleek()   +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title=element_blank(),
        legend.position="none") +
  scale_x_continuous(limits = c(0, 1))


par.long <- par.posts |>
  mutate(CU_f = CU,
         alpha = exp(ln_a)) |>
  left_join(CU_name_lookup, by="CU_f")

a <- ggplot(par.long, aes(alpha, fill = CU_pretty, color = CU_pretty)) +
  geom_density(alpha = 0.3,bw=0.4) +
  theme(legend.position = "bottom") +
  labs(x = "Intrinsic productivity", y = "Posterior density") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_sleek()   +
  theme(legend.position = c(0.8,0.625),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title=element_blank()) +
  scale_x_continuous(limits = c(0, 15))

cowplot::plot_grid(a, b, c, labels="auto", ncol=1)

ggsave(here("csasdown/figure/par-ref-hist.PNG"), width = 675*2, height = 975*2, 
       units="px", dpi=240)


ggplot(par.long, aes(log(alpha), fill = CU_pretty, color = CU_pretty)) +
  geom_density(alpha = 0.3,bw=0.4) +
  theme(legend.position = "bottom") +
  labs(x = "Intrinsic productivity", y = "Posterior density") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_sleek()   +
  theme(legend.position = c(0.8,0.625),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title=element_blank()) + 
  scale_x_continuous(limits = c(-2, 6))
