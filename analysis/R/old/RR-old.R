# scripts to fit multi-CU run-reconstruction model to GSI and border passage observations, then wrangle output to create files for SR models

library(here)
library(tidyverse)
source(here("analysis/R/run-reconstructions/initRR.R"))

# generate data file if inputs updated
# source(here("analysis/R/run-reconstructions/procData.R"))
# processData()

# fit run-reconstruction model
rpt <- fitRR()

# process RR model output and save time-series of spawner abundance
cdn_harvest <- read.csv(here("analysis/data/raw/YkCk_Harvest_CA_Data_2024.csv")) %>%
  filter(Type == "CA_Mainstem",
         Year > 1984) %>%
  select(Year, Estimate)

border_passage <- rpt[["runSize_t"]]
cdn_er <- cdn_harvest$Estimate/(cdn_harvest$Estimate+border_passage)

CU_border_passage <- exp(rpt$lnRunSize_st)*1e-3
CU_spwn <- CU_border_passage * (1-cdn_er)

colnames(CU_spwn) <- seq(1985,2024)
rownames(CU_spwn) <- c("NorthernYukonR.andtribs.","Whiteandtribs.","Pelly","Stewart","Nordenskiold","YukonR.Teslinheadwaters","MiddleYukonR.andtribs.","UpperYukonR.")

# estimate proportional contribution of Big Salmon to Middle CU, generate full time series, remove from Middle CU
big_salmon <- read.csv(here("analysis/data/raw/trib-spwn.csv")) |>
  filter(system == "bigsalmon",
         year < 2025)

prop_comp <- mean(big_salmon$estimate/(CU_spwn[7,c(21:40)]*1000), na.rm=TRUE)

big_salmon_inf <- CU_spwn[7,]*1000 * prop_comp
big_salmon_inf_2 <- big_salmon_inf
big_salmon_inf_2[c(21:40)] <- big_salmon$estimate
big_salmon_inf_2[38] <- big_salmon_inf[38]
big_salmon_recon <- big_salmon_inf_2/1000
CU_spwn[7,] <- CU_spwn[7,]-big_salmon_recon

# turn spawners into long data frame
CU_spwn_df <- as.data.frame(CU_spwn)
CU_spwn_df$CU <- rownames(CU_spwn)

CU_spawn_long <- CU_spwn_df %>% 
  pivot_longer(!CU, names_to = 'Year', values_to = 'spawn')

Rse <- filter(rpt$sdrpt,par=="runSize_st")
CU_spawn_long$lwr <- Rse$lCI*1e-3* (1-cdn_er)
CU_spawn_long$upr <- Rse$uCI*1e-3* (1-cdn_er)
CU_spawn_long$se <- Rse$se*1e-3

mssr_spwn <- CU_spawn_long %>%
  mutate(year = as.numeric(Year),
         cv = se/spawn,
         stock = CU,
         mean = as.numeric(spawn)*1000,
         obs = 1) %>%
  select(stock, year, mean, se, cv, obs,lwr, upr)

# add Big Salmon into spawner dataframe
year <- seq(1985,2024)
big.S <- cbind(year,(big_salmon_recon*1000))
big.S <- as.data.frame(big.S)
big.S$stock <- "Big.Salmon"
big.S$se <- "na"
big.S$cv <- 0.50
big.S$obs <- 1
big.S$cv[c(21:40)] <- 0.05
big.S <- big.S[,c(3,1,2,4,5,6)]
big.S$lwr <- "na"
big.S$upr <- "na"

colnames(big.S)[3] <- "mean"

mssr_spwn_2 <- rbind(mssr_spwn,big.S)

write.csv(mssr_spwn_2, here("analysis/data/raw/esc-data.csv"),row.names = F)

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
write.csv(harvest,here("analysis/data/raw/harvest-data.csv"))
