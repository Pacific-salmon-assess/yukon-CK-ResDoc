library(tidyverse)
library(here)
library(ggsidekick)

# harvest by AK district and year
AK_dist_harv <- read.csv(here("analysis/data/raw/YkCk_Harvest_byDistrTypeStockAge.csv"))

# CDN harvest by district and year
CDN_AK_dist_harv <- AK_dist_harv |>
  filter(Stock == "Upper",
         District != 7) |>
  group_by(Year, District, Fishery) |>
  summarise(total_harvest = sum(Total)) |>
  filter(Year %in% c(2005,2006,2007,2008,2009,2010))


# harvest genetics by AK district and year
AK_dist_harv_genetics <- read.csv(here("analysis/data/raw/AK-harvest-genetics.csv")) |>
  group_by(Year, District, Fishery, Stock) |>
  summarise(avg_contr = mean(Contribution))|>
  na.omit()

AK_dist_harv_genetics_total_contr <- read.csv(here("analysis/data/raw/AK-harvest-genetics.csv")) |>
  group_by(Year, District, Fishery, Stock) |>
  summarise(avg_contr = mean(Contribution)) |>
  group_by(Year, District, Fishery) |>
  summarise(total_contr = sum(avg_contr)) |>
  na.omit()

AK_harv_genetics <- AK_dist_harv_genetics |>
  left_join(AK_dist_harv_genetics_total_contr, by=c("Year", "District", "Fishery")) |>
  mutate(contr = avg_contr/total_contr) |>
  select(Year, District, Fishery, Stock, contr) |>
  left_join(CDN_AK_dist_harv, by=c("Year", "District", "Fishery")) |>
  mutate(harvest = contr*total_harvest) |>
  group_by(Year, Stock) |>
  summarize(stock_harvest = sum(harvest)) |>
  mutate(harv_per_contr = stock_harvest/sum(stock_harvest))

# border comps by CU and year
CU_border_comps <- read.csv(here("analysis/data/generated/CU-border-comps.csv")) |>
  mutate(Stock = case_when(
           CU_f == "NorthernYukonR.andtribs." ~ "Border",
           CU_f == "Whiteandtribs." ~ "Border",
           CU_f == "Pelly" ~ "Pelly",
           CU_f == "Stewart" ~ "Pelly",
           CU_f == "Nordenskiold" ~ "Carmacks",
           CU_f == "MiddleYukonR.andtribs." ~ "Carmacks",
           CU_f == "YukonR.Teslinheadwaters" ~ "Carmacks",
           CU_f == "Big.Salmon" ~ "Carmacks",
           CU_f == "UpperYukonR." ~ "Takhini")) |>
  rename(Year=year) |>
  group_by(Year,Stock) |>
  summarize(RR_per_contr=(sum(contr)/100)) |>
  filter(Year %in% c(2005,2006,2007,2008,2009,2010))

fish_select_explore <- AK_harv_genetics |>
  left_join(CU_border_comps, by=c("Year","Stock")) |>
  select(!stock_harvest) |>
  rename(Harvest = harv_per_contr) |>
  rename(Border = RR_per_contr) |>
  pivot_longer(!c(Year,Stock), names_to = "Source", values_to = "Contribution" )

ggplot(fish_select_explore, aes(x = Stock, y = Contribution, fill = Source)) +
  geom_bar(position="dodge", stat = "identity") +
  theme_sleek() +
  labs(x = "CU group", y = "Contribution") +
  facet_wrap(~Year)

my.ggsave(here("analysis/plots/RR/harbest_vs_border_composition.PNG"), height=4.25, width=8)
my.ggsave(here("csasdown/figure/harbest_vs_border_composition.PNG"), height=4.25, width=8)

