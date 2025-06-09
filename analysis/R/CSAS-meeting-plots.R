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

# compare alt disperson RR models ----

esc_alts <- read.csv(here("analysis/data/generated/esc-data_alt_disp.csv")) 

ggplot(esc_alts,aes(x = year, y = mean, col=as.factor(dispersion))) + 
  geom_line(lwd = 0.8) +
  xlab("Year") +
  ylab("spawners (000s)")+
  theme_sleek() +
  theme(strip.text = element_text(size=10)) +
  facet_wrap(~cu, scales = "free_y")

my.ggsave(here("analysis/plots/RR/cu-esc-alt-RR-disp.PNG"))

# compare CU and agg Seq ---

posteriors <- read.csv(here("analysis/data/generated/AR1_posteriors.csv"))
posteriors$seq <- posteriors$ln_a/posteriors$beta

X<-matrix(NA,nrow=1000,ncol=1)
for(i in 1:1000){
x  <- posteriors|>
  group_by(CU) |>
  slice_sample() |>
  ungroup() |>
  summarize(agg_seq = sum(seq))
X[i] <- x$agg_seq
}

agg_seq <- as.data.frame(X); colnames(agg_seq) <- "seq"
ggplot(agg_seq, aes(x = seq)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density()

summarise(agg_seq)

median(agg_seq$seq)

# compare ref points across alt RR models ----


# read in data ---------------------------------------------------------------------------
# model fits --- 
AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1(disp=0.02)"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

# process data and fits to make plots later ----------------------------------------------
bench.posts <- NULL
bench.posts.all <- NULL
for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs 
  
  # AR1 (base S-R) models -----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  bench <- matrix(NA,length(sub_pars$lnalpha),6,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("lnalpha","Sgen", "Smsy", "Umsy", "Seq", "Smsr")))
  
  # get benchmarks & pars (AR1 model)------------------------------------------------------------------
  
  for(j in 1:length(sub_pars$lnalpha)){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    bench[j,1] <- ln_a
    bench[j,3] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,2] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,3]) #S_gen
    bench[j,4] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,5] <- ln_a/b #S_eq
    bench[j,6] <- 1/b #S_msr 
  }
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))

}  # End data wrangling loop by CU 

bench.posts$disp <- 0.02
bench.posts.all <- rbind(bench.posts.all, as.data.frame(bench.posts))
      
AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1(disp=0.1)"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

# process data and fits to make plots later ----------------------------------------------
bench.posts <- NULL

for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs 
  
  # AR1 (base S-R) models -----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  bench <- matrix(NA,length(sub_pars$lnalpha),6,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("lnalpha","Sgen", "Smsy", "Umsy", "Seq", "Smsr")))
  
  # get benchmarks & pars (AR1 model)------------------------------------------------------------------
  
  for(j in 1:length(sub_pars$lnalpha)){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    bench[j,1] <- ln_a
    bench[j,3] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,2] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,3]) #S_gen
    bench[j,4] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,5] <- ln_a/b #S_eq
    bench[j,6] <- 1/b #S_msr 
  }
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
}  # End data wrangling loop by CU 

bench.posts$disp <- 0.1
bench.posts.all <- rbind(bench.posts.all, as.data.frame(bench.posts))

AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1(disp=0.5)"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

# process data and fits to make plots later ----------------------------------------------
bench.posts <- NULL

for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs 
  
  # AR1 (base S-R) models -----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  bench <- matrix(NA,length(sub_pars$lnalpha),6,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("lnalpha","Sgen", "Smsy", "Umsy", "Seq", "Smsr")))
  
  # get benchmarks & pars (AR1 model)------------------------------------------------------------------
  
  for(j in 1:length(sub_pars$lnalpha)){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    bench[j,1] <- ln_a
    bench[j,3] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,2] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,3]) #S_gen
    bench[j,4] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,5] <- ln_a/b #S_eq
    bench[j,6] <- 1/b #S_msr 
  }
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
}  # End data wrangling loop by CU 

bench.posts$disp <- 0.5
bench.posts.all <- rbind(bench.posts.all, as.data.frame(bench.posts))

AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1(disp=1)"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

# process data and fits to make plots later ----------------------------------------------
bench.posts <- NULL

for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs 
  
  # AR1 (base S-R) models -----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  bench <- matrix(NA,length(sub_pars$lnalpha),6,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("lnalpha","Sgen", "Smsy", "Umsy", "Seq", "Smsr")))
  
  # get benchmarks & pars (AR1 model)------------------------------------------------------------------
  
  for(j in 1:length(sub_pars$lnalpha)){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    bench[j,1] <- ln_a
    bench[j,3] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,2] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,3]) #S_gen
    bench[j,4] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,5] <- ln_a/b #S_eq
    bench[j,6] <- 1/b #S_msr 
  }
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
}  # End data wrangling loop by CU 

bench.posts$disp <- 1
bench.posts.all <- rbind(bench.posts.all, as.data.frame(bench.posts))

AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1(disp=1.5)"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

# process data and fits to make plots later ----------------------------------------------
bench.posts <- NULL

for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs 
  
  # AR1 (base S-R) models -----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  bench <- matrix(NA,length(sub_pars$lnalpha),6,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("lnalpha","Sgen", "Smsy", "Umsy", "Seq", "Smsr")))
  
  # get benchmarks & pars (AR1 model)------------------------------------------------------------------
  
  for(j in 1:length(sub_pars$lnalpha)){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    bench[j,1] <- ln_a
    bench[j,3] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,2] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,3]) #S_gen
    bench[j,4] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,5] <- ln_a/b #S_eq
    bench[j,6] <- 1/b #S_msr 
  }
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
}  # End data wrangling loop by CU 

bench.posts$disp <- 1.5
bench.posts.all <- rbind(bench.posts.all, as.data.frame(bench.posts))


AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1(disp=2)"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

# process data and fits to make plots later ----------------------------------------------
bench.posts <- NULL

for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs 
  
  # AR1 (base S-R) models -----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  bench <- matrix(NA,length(sub_pars$lnalpha),6,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("lnalpha","Sgen", "Smsy", "Umsy", "Seq", "Smsr")))
  
  # get benchmarks & pars (AR1 model)------------------------------------------------------------------
  
  for(j in 1:length(sub_pars$lnalpha)){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    bench[j,1] <- ln_a
    bench[j,3] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,2] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,3]) #S_gen
    bench[j,4] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,5] <- ln_a/b #S_eq
    bench[j,6] <- 1/b #S_msr 
  }
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
}  # End data wrangling loop by CU 

bench.posts$disp <- 2
bench.posts.all <- rbind(bench.posts.all, as.data.frame(bench.posts))

write.csv(bench.posts.all, here("analysis/data/generated/bench.posts.all.csv"), 
          row.names = FALSE)

a <- ggplot(bench.posts.all, aes(x = CU, y = lnalpha, fill = as.factor(disp))) +
  geom_boxplot(outlier.shape=NA) +
  theme_sleek() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank()) +
  labs(x = "Conservation Unit", y = expression(Log(alpha))) +
  ylim(0,4)

b <- ggplot(bench.posts.all, aes(x = CU, y = Smsr/1000, fill = as.factor(disp))) +
  geom_boxplot(outlier.shape=NA) +
  theme_sleek() +
  theme(legend.position = c(0.9,0.625),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank()) +
  labs(x = "Conservation Unit", y = bquote(S[MSR]~(`000s`))) +
  ylim(0,50)



c <- ggplot(bench.posts.all, aes(x = CU, y = Seq/1000, fill = as.factor(disp))) +
  geom_boxplot(outlier.shape=NA) +
  theme_sleek() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank()) +
  labs(x = "Conservation Unit", y = bquote(S[EQ]~(`000s`))) +
  ylim(0,50)

cowplot::plot_grid(a, b, c, labels="auto", ncol=1,rel_heights=c(0.6,0.6,1))

ggsave(here("analysis/plots/SR-models/par-compare-alt-RR.PNG"), width=700*2, height=900*2, dpi=240,
       units="px")


# summarize ----

bench.posts.all |>
  pivot_longer(cols=lnalpha:Smsr,names_to="parameter", values_to="value") |>
  group_by(parameter, CU) |>
  summarize(avg=mean(value))
