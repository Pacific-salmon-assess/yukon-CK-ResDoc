# scripts to fit multi-CU run-reconstruction model to GSI and border passage observations, then wrangle output to create files for SR models

library(here)
library(tidyverse)
source(here("analysis/R/run-reconstructions/initRR.R"))

# generate data file if inputs updated
# source(here("analysis/R/run-reconstructions/procData.R"))
# processData()

# fit run-reconstruction model
# rpt <- fitRR()

# load run-reconstruction model fit

load(here("analysis/R/run-reconstructions/fittedMod/rpt.Rdata"))

# process RR model output and save time-series of spawner abundance
cdn_harvest <- read.csv(here("analysis/data/raw/YkCk_Harvest_CA_Data_2024.csv")) %>%
  filter(Type == "CA_Mainstem",
         Year > 1984) %>%
  select(Year, Estimate)

border_passage <- rpt[["runSize_t"]]
cdn_er <- cdn_harvest$Estimate/(cdn_harvest$Estimate+border_passage)

CU_border_passage <- exp(rpt$lnRunSize_st)
CU_cdn_er <- matrix( data=rep((1-cdn_er),each=8), nrow=8, ncol=40)
CU_spwn <- CU_border_passage * CU_cdn_er

# extract CIs for border passage, calculate spawners
Rse <- filter(rpt$sdrpt,par=="runSize_st")
CU_spwn_upr <- matrix( data=Rse$uCI, nrow=8, ncol=40) * CU_cdn_er
CU_spwn_lwr <- matrix( data=Rse$lCI, nrow=8, ncol=40) * CU_cdn_er

# add CUs and years
colnames(CU_spwn) <- seq(1985,2024)
rownames(CU_spwn) <- c("NorthernYukonR.andtribs.","Whiteandtribs.","Pelly","Stewart","Nordenskiold","YukonR.Teslinheadwaters","MiddleYukonR.andtribs.","UpperYukonR.")

colnames(CU_spwn_upr) <- seq(1985,2024)
rownames(CU_spwn_upr) <- c("NorthernYukonR.andtribs.","Whiteandtribs.","Pelly","Stewart","Nordenskiold","YukonR.Teslinheadwaters","MiddleYukonR.andtribs.","UpperYukonR.")

colnames(CU_spwn_lwr) <- seq(1985,2024)
rownames(CU_spwn_lwr) <- c("NorthernYukonR.andtribs.","Whiteandtribs.","Pelly","Stewart","Nordenskiold","YukonR.Teslinheadwaters","MiddleYukonR.andtribs.","UpperYukonR.")

# estimate proportional contribution of Big Salmon to Middle CU, generate full time series, remove from Middle CU
big_salmon <- read.csv(here("analysis/data/raw/trib-spwn.csv")) |>
  filter(system == "bigsalmon",
         type == "sonar",
         year < 2025)

prop_comp <- mean(big_salmon$estimate/(CU_spwn[7,c(21:40)]), na.rm=TRUE)
prop_comp_ind <- big_salmon$estimate/CU_spwn[7,c(21:40)]
hist(prop_comp_ind)
mu_bs_prop <- mean(prop_comp_ind, na.rm=TRUE)
sd_bs_prop <- sd(prop_comp_ind, na.rm=TRUE)

r_bs_prop <- rnorm(1000,mu_bs_prop,sd_bs_prop)
r_bs_prop <- abs(r_bs_prop)

mu_my <- CU_spwn[7,]
sd_my <- (CU_spwn_upr[7,]-CU_spwn_lwr[7,])/4

bs <- matrix(NA,1000,40)
my <- matrix(NA,1000,40)

for (i in 1:1000){
  r_bs_prop <- rnorm(1,mu_bs_prop,sd_bs_prop)
  r_my <- rnorm(mu_my,mean=mu_my,sd=sd_my)
  bs[i,] <- r_my*r_bs_prop
  my[i,] <- r_my-bs[i,]
}  
  
bs_mean <- apply(bs,2,mean)
bs_upr <- apply(bs,2,quantile, probs=c(0.95))
bs_lwr <- apply(bs,2,quantile, probs=c(0.05))
bs_cv <- apply(bs,2,sd)/bs_mean

my_mean <- apply(my,2,mean)
my_cv <- apply(my,2,sd)/my_mean
my_upr <- apply(my,2,quantile, probs=c(0.95))
my_lwr <- apply(my,2,quantile, probs=c(0.05))

big_salmon_inf_2 <- bs_mean
big_salmon_inf_2[c(21:40)] <- big_salmon$estimate
big_salmon_inf_2[38] <- bs_mean[38]
big_salmon_recon <- big_salmon_inf_2
CU_spwn[7,] <- my_mean
CU_spwn_upr[7,] <- my_upr
CU_spwn_lwr[7,] <- my_lwr

# estimate hatchery contribution from Whitehorse fishway, remove from Upper Yukon CU
hatch_upper <- read.csv(here("analysis/data/raw/trib-spwn.csv")) |>
  filter(system == "whitehorse",
         year > 1984) |>
  mutate(hatch = estimate*hatch_contrib)

CU_spwn[8,] <- CU_spwn[8,]-hatch_upper$hatch
CU_spwn_upr[8,] <- CU_spwn_upr[8,]-hatch_upper$hatch
CU_spwn_lwr[8,] <- CU_spwn_lwr[8,]-hatch_upper$hatch

# turn spawners into long data frame
CU_spwn_df <- as.data.frame(CU_spwn)
CU_spwn_df$CU <- rownames(CU_spwn)

CU_spawn_long <- CU_spwn_df %>% 
  pivot_longer(!CU, names_to = 'Year', values_to = 'spawn')

CU_spwn_lwr_df <- as.data.frame(CU_spwn_lwr)
CU_spwn_lwr_df$CU <- rownames(CU_spwn_lwr)

CU_spawn_lwr_long <- CU_spwn_lwr_df %>% 
  pivot_longer(!CU, names_to = 'Year', values_to = 'spawn.lwr') |>
  mutate(spawn.lwr = ifelse(spawn.lwr < 0, 10, spawn.lwr))

CU_spwn_upr_df <- as.data.frame(CU_spwn_upr)
CU_spwn_upr_df$CU <- rownames(CU_spwn_upr)

CU_spawn_upr_long <- CU_spwn_upr_df %>% 
  pivot_longer(!CU, names_to = 'Year', values_to = 'spawn.upr') 

CU_spawn_2 <- cbind(CU_spawn_long,CU_spawn_lwr_long[,3],CU_spawn_upr_long[,3]) |>
  mutate(cv = ((spawn.upr-spawn.lwr)/4)/spawn) |>
  rename(stock = CU, year = Year, mean = spawn, upper = spawn.upr, lower = spawn.lwr)

# add Big Salmon into spawner dataframe
year <- seq(1985,2024)
big.S <- cbind(year,(big_salmon_recon))
big.S <- as.data.frame(big.S)
big.S$stock <- "Big.Salmon"
big.S$cv <- bs_cv
big.S$cv[c(21:40)] <- 0.05
big.S$lwr <- bs_lwr
big.S$upr <- bs_upr

big.S$lwr[21:40] <- (big_salmon_recon+(big_salmon_recon*0.1))[21:40]
big.S$upr[21:40]  <- (big_salmon_recon-(big_salmon_recon*0.1))[c(21:40)]
big.S$lwr[38] <- bs_lwr[38]
big.S$upr[38] <- bs_upr[38]

big.S <- big.S |>
  rename(mean = V2) |>
  mutate(lower = lwr,
         upper = upr) |>
  select(stock, year, mean, upper, lower, cv) 

# bind data frames
mssr_spwn_2 <- rbind(CU_spawn_2,big.S)

write.csv(mssr_spwn_2, here("analysis/data/generated/esc-data.csv"),row.names = F)

# calculate CU specific harvest based on reconstructed spawner abundance and aggregate exploitation rate
er <- read.csv(here("analysis/data/raw/rr-table.csv"))

cdn_er <- er %>%
  filter(stock == "Canadian",
         Year > 1984) %>%
  mutate(er = Harvest.rate..../100,
         cv = Harvest.CV)

agg_er <- rep(as.vector(cdn_er$er),9)
agg_cv <- rep(as.vector(cdn_er$cv),9)

harv <- cbind(mssr_spwn_2, agg_er,agg_cv)
harv$harv <- (harv$mean/(1-harv$agg_er))*harv$agg_er
harvest <- harv[,c(1,2,9,8)]
colnames(harvest) <- c("population", "year","harv","cv")
write.csv(harvest,here("analysis/data/generated/harvest-data.csv"))






















