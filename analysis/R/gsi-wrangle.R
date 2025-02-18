# wrangle raw GSI into long format for RR model

library(tidyverse)
library(here)

raw.gsi <- read.csv(here("analysis/data/raw/Final_Combined_results_2025-02-17_ER_AT_JS.csv"))

long.gsi <- raw.gsi %>% pivot_longer(
    cols = NYR:UpperYR,
    names_to = "CU",
    values_to = "prob") %>%
  filter(use == "YES") %>%
  rename(julian = Day,
         year = Year,
         fish = Fish) %>%
  mutate(CU_no = case_when(
    CU == "NYR" ~ 76,
    CU == "White" ~ 75,
    CU == "Pelly" ~ 72,
    CU == "Stew" ~ 74,
    CU == "Norden" ~ 71,
    CU == "Teslin" ~ 68,
    CU == "MidYR" ~ 73,
    CU == "UpperYR" ~ 69)) %>%
  select(year, gear, julian, fish, CU_no, CU, prob)

write.csv(long.gsi,"analysis/data/raew/border-gsi-table-2024-update.csv")
