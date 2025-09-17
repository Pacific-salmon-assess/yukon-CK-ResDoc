library(here)
library(tidyverse)

esc <- read.csv(here("analysis/data/generated/esc-data.csv")) |>
  mutate_at(2:6, as.numeric)

spwn_cosewic <- esc |>
  mutate(CUID = case_when(stock == "NorthernYukonR.andtribs." ~ "1207",
                          stock == "Whiteandtribs." ~ "1206",
                          stock == "Pelly" ~ "1203",
                          stock == "Stewart" ~ "1205",
                          stock == "Big.Salmon" ~ "1201",
                          stock == "YukonR.Teslinheadwaters" ~ "1212",
                          stock == "Nordenskiold" ~ "1202",
                          stock == "MiddleYukonR.andtribs." ~ "1204",
                          stock == "UpperYukonR." ~ "1211"),
         Year = year) %>%
  group_by(CUID, Year)%>%
  summarise(Spawner.Abundance = mean(mean)) %>%
  mutate(Species = "Chinook",
         Region = "Yukon") %>%
  select(CUID, Species, Region, Year, Spawner.Abundance)

write.csv(spwn_cosewic, here("analysis/data/generated/spwn.cosewic.04Sep2025.csv"),row.names = F)
