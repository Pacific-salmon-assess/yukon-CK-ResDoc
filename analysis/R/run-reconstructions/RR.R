# scripts to fit multi-CU run-reconstruction model to GSI and border passage observations, then wrangle output to create files for SR models

library(here)
library(tidyverse)
source(here("analysis/R/run-reconstructions/initRR.R"))

# generate data file if inputs updated
# source(here("analysis/R/run-reconstructions/procData.R"))
# processData()

# fit run-reconstruction model
#rpt <- fitRR()

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

big_salmon_inf <- CU_spwn[7,] * prop_comp
big_salmon_inf_2 <- big_salmon_inf
big_salmon_inf_2[c(21:40)] <- big_salmon$estimate
big_salmon_inf_2[38] <- big_salmon_inf[38]
big_salmon_recon <- big_salmon_inf_2
CU_spwn[7,] <- CU_spwn[7,]-big_salmon_recon
CU_spwn_upr[7,] <- CU_spwn_upr[7,]-big_salmon_recon
CU_spwn_lwr[7,] <- CU_spwn_lwr[7,]-big_salmon_recon

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
big.S$cv <- 0.25
big.S$cv[c(21:40)] <- 0.05
big.S$lwr = "na"
big.S$upr = "na"

big.S <- big.S |>
  rename(mean = V2) |>
  mutate(lower = (mean - (2*(mean*cv))),
         upper = (mean + (2*(mean*cv)))) |>
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

























