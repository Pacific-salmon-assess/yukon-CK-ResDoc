# histograms of reference points ----
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

# Comapre benchmarks ----
ref_points <- read.csv(here("analysis/data/generated/bench_par_table.csv")) 
dem_ref_points <- read.csv(here("analysis/data/generated/demographic_parameters.csv")) 

simple <- ref_points |>
  filter(bench.par %in% c("Sgen","Smsr.20", "Smsr.40","Smsy","Seq","Smsr")) |>
  mutate(median = X50.,
         lwr=X10.,
         upr=X90.,
         par=bench.par,
         model="spawners",
         mecase_when()) |>
  select(CU,par,median,lwr,upr,model)

simple.1 <- simple |>
  filter(par == "Smsy") |>
  mutate(par = "Smsy.80",
         median = median*0.8,
         lwr=lwr*0.8,
         upr=upr*0.8) |>
  select(CU,par,median,lwr,upr,model)

demo.1 <- dem_ref_points |>
  filter(period == "recent",
         par == "Smsr") |>
  mutate(par = "Smsr.egg-mass",
         model="egg-mass",
         lwr=lower,
         upr=upper) |>
  select(CU,par,median,lwr,upr,model)

demo.2 <- demo.1 |>
  mutate(par = "Smsr.20.egg-mass",
         median = median*0.2,
         lwr=lwr*0.2,
         upr=upr*0.2) |>
  select(CU,par,median,lwr,upr,model)

demo.3 <- demo.1 |>
  mutate(par = "Smsr.40.egg-mass",
         median = median*0.4,
         lwr=lwr*0.4,
         upr=upr*0.4) |>
  select(CU,par,median,lwr,upr,model)

demo <- rbind(demo.1,demo.2,demo.3)

pars <- rbind(simple,simple.1,demo) 

a <- ggplot(pars |>filter(
  par %in% c("Sgen","Smsr.20","Smsr.20.egg-mass")), aes(x = CU, y = median, fill = par)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(aes(ymin = lwr, ymax = upr,col = par), width = 0,position=position_dodge(0.9)) +
  theme_sleek() +
  theme(legend.position = c(0.8,0.825),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank()) +
  labs(x = "Conservation Unit", y = "Spawners") 


pars2 <- pars |>filter(
  par %in% c("Smsy.80","Smsr.40","Smsr.40.egg-mass"))
pars2$par2 <- factor(pars2$par,levels=c("Smsy.80","Smsr.40","Smsr.40.egg-mass"))

b <- ggplot(pars2, aes(x = CU, y = median, fill = par2)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(aes(ymin = lwr, ymax = upr,col = par2), width = 0,position=position_dodge(0.9)) +
  theme_sleek()  +
  theme(legend.position = c(0.8,0.825),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank()) +
  labs(x = "Conservation Unit", y = "Spawners") 

c <- ggplot(pars |>filter(
  par %in% c("Seq","Smsr","Smsr.egg-mass")), aes(x = CU, y = median, fill = par)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(aes(ymin = lwr, ymax = upr,col = par), width = 0,position=position_dodge(0.9)) +
  theme_sleek() +
  theme(legend.position = c(0.8,0.825),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank()) +
  labs(x = "Conservation Unit", y = "Spawners") 

cowplot::plot_grid(a, b, c, labels="auto", ncol=1,rel_heights=c(0.6,0.6,1))

ggsave(here("analysis/plots/SR-models/bench-compare.PNG"), width=700*2, height=900*2, dpi=240,
       units="px")
