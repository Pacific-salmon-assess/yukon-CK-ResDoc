#wrangling ASL data i

library(tidyverse)
library(here)
library(zoo)

# Hamachan selectivity corrected age and sex comps 1981-2006
fw_age_sex <- read.csv(here("analysis/data/raw/fw-border-age-sex-comps.csv"))

# Hamachan length raw data with selectivity 1981-2006
fw_l <- read.csv(here("analysis/data/raw/FWborder.csv"))

# Eage ASL data 2005-2020
eagle_age_sex_len <- read.csv(here("analysis/data/raw/ASL_Output_Chinook_Eagle_2005-2024.csv"))

## wrangling for full age and sex dataset
#  eagle ASL prep
eagle_asl <- eagle_age_sex_len %>% 
  filter(!is.na(Total.Age), Sex != "unknown") |>
  group_by(Sample.Year, Sex, Total.Age) |>
  count() |>
  group_by(Sample.Year) |>
  mutate(total_n=sum(n)) |>
  ungroup() |>
  mutate(prop_n=n/total_n)

eagle_as_2005_2024 <- eagle_asl |>
  select(Sample.Year,Sex, Total.Age, prop_n) |>
  pivot_wider(names_from = Total.Age,
              values_from =prop_n) |>
  rename(age_3="3",
         age_4="4",
         age_5="5",
         age_6="6",
         age_7="7",
         age_8="8")%>%
  mutate_if(is.numeric , replace_na, replace = 0) |>
  mutate(age_4=age_3+age_4,
         age_7=age_7+age_8) |>
  filter(Sample.Year>2006) |>
  select(Sample.Year,Sex, age_4, age_5, age_6, age_7)|>
  filter(Sample.Year<2025)

#fish wheel ASL prep
fw_as_1984_2006 <- fw_age_sex |>
  mutate_if(is.numeric , replace_na, replace = 0) |>
  rename(Sample.Year=Year,
         Sex=sex) |>
  mutate(age_4=X3_yr+X4_yr,
         age_5=X5_yr,
         age_6=X6_yr,
         age_7=X7_yr+X8_yr) |>
  filter(Sample.Year>1984) |>
    select(Sample.Year,Sex, age_4, age_5, age_6, age_7)

# merge datasets
cdn_as_1985_2024 <- rbind(fw_as_1984_2006,eagle_as_2005_2024)

write.csv(cdn_as_1985_2024[which(cdn_as_1985_2024$Sex=="female"),], here("analysis/data/raw/female_age_comps.csv"),row.names = FALSE)

## wrangling for full length dataset
# fish wheel data
fw_fl <- fw_l |> 
  filter(!is.na(Age)) |>
  group_by(Year, Sex, Age) |>
  summarise(mean_fl=sum(Length*nsel)/sum(nsel)) |>
  pivot_wider(names_from = Age,
              values_from = mean_fl) |>
  rename(age_4="4",
         age_5="5",
         age_6="6",
         age_7="7",
         Sample.Year=Year) |>
  select(Sample.Year,Sex, age_4, age_5, age_6, age_7) |>
  filter(!is.na(Sex),Sex!= "male",Sample.Year<2007) |>
  ungroup()

# fill years and ages with no data with next year observation
fw_fl$age_4 <- na.locf(zoo(fw_fl$age_4), fromLast = TRUE)
fw_fl$age_5 <- na.locf(zoo(fw_fl$age_5), fromLast = TRUE)
fw_fl$age_6 <- na.locf(zoo(fw_fl$age_6), fromLast = TRUE)
fw_fl$age_7 <- na.locf(zoo(fw_fl$age_7), fromLast = TRUE)

#fw_fl <- fw_fl[,-1]

# eagle data
eagle_fl <- eagle_age_sex_len |> 
  filter(!is.na(Total.Age), Sex != "unknown") |>
  group_by(Sample.Year, Sex, Total.Age) |>
  summarise(fl=mean(Length, na.rm = TRUE)) |>
  select(Sample.Year, Sex, Total.Age, fl) |>
  pivot_wider(names_from = Total.Age,
              values_from = fl) |>
  rename(age_4="4",
         age_5="5",
         age_6="6",
         age_7="7") |>
         select(Sample.Year, Sex, age_4, age_5, age_6, age_7) |>
  filter(Sex!= "male", Sample.Year>2006) |>
  ungroup()

eagle_fl$age_7 <- na.locf(zoo(eagle_fl$age_7), fromfirst = TRUE)

# merge datasets, add in "data" for missing 1984 year
cdn_len_1985_2024 <- rbind(fw_fl,eagle_fl) |>
  filter(Sample.Year > 1984)

write.csv(cdn_len_1985_2024, here("analysis/data/raw/female_length_comps.csv"),row.names = FALSE)

#convert to expected eggs, or egg mass, per female in each age class assuming average length
female_length_comps_eggs <- 9.35e-4*cdn_len_1985_2024[,3:6]^2.36
female_length_comps_egg_mass <- 8.71e-12*cdn_len_1985_2024[,3:6]^4.83

write.csv(female_length_comps_eggs, here("analysis/data/raw/female_length_comps_eggs.csv"),row.names = FALSE)
write.csv(female_length_comps_egg_mass, here("analysis/data/raw/female_length_comps_egg_mass.csv"),row.names = FALSE)

