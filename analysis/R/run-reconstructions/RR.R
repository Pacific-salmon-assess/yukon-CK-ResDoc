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
cdn_harvest <- read.csv(here("analysis/data/raw/YkCk_Harvest_CA_Data.csv")) %>%
  filter(Type == "CA_Mainstem",
         Year > 1984) %>%
  select(Year, Estimate)

border_passage <- rpt[["runSize_t"]]
cdn_er <- cdn_harvest$Estimate/(cdn_harvest$Estimate+border_passage)

CU_border_passage <- exp(rpt$lnRunSize_st)*1e-3
CU_spwn <- CU_border_passage * (1-cdn_er)

colnames(CU_spwn) <- seq(1985,2023)
rownames(CU_spwn) <- c("NorthernYukonR.andtribs.","Whiteandtribs.","Pelly","Stewart","Nordenskiold","YukonR.Teslinheadwaters","MiddleYukonR.andtribs.","UpperYukonR.")

CU_spwn_df <- as.data.frame(CU_spwn)
CU_spwn_df$CU <- rownames(CU_spwn)

CU_spawn_long <- CU_spwn_df %>% 
  pivot_longer(!CU, names_to = 'Year', values_to = 'spawn')

Rse <- filter(rpt$sdrpt,par=="runSize_st")
CU_spawn_long$lwr <- Rse$lCI*1e-3
CU_spawn_long$upr <- Rse$uCI*1e-3
CU_spawn_long$se <- Rse$se*1e-3

mssr_spwn <- CU_spawn_long %>%
  mutate(year = as.numeric(Year),
         cv = se/spawn,
         stock = CU,
         mean = as.numeric(spawn)*1000,
         obs = 1) %>%
  select(stock, year, mean, se, cv, obs)

write.csv(mssr_spwn, here("analysis/data/raw/esc-data.csv"),row.names = F)

# calculate CU specific harvest based on reconstructed spawner abundance and aggregate exploitation rate
er <- read.csv(here("analysis/data/raw/rr-table.csv"))

cdn_er <- er %>%
  filter(stock == "Canadian",
         Year > 1984) %>%
  mutate(er = Harvest.rate..../100,
         cv = Harvest.CV)

agg_er <- rep(as.vector(cdn_er$er),8)
agg_cv <- rep(as.vector(cdn_er$cv),8)

harv <- cbind(mssr_spwn, agg_er,agg_cv)
harv$harv <- (harv$mean/(1-harv$agg_er))*harv$agg_er
harvest <- harv[,c(1,2,9,8)]
colnames(harvest) <- c("population", "year","harv","cv")
write.csv(harvest,here("analysis/data/raw/harvest-data.csv"))
