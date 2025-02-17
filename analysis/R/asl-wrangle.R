#wrangling ASL data into 2 matrices

library(tidyverse)
library(here)

#helper funs-----------------------------------------------------------------------------
#needed since using native pipe (i.e. referencing '.' doesn't work)
my_replace <- function(x){
  return(replace(
    x = x,
    list = is.na(x),
    values = 0))
}

#read in and wrangle asl-------------------------------------------------------------- 
asl <- read.csv(here("analysis/data/raw/ASL_Output_Chinook_Eagle_2005-2019.csv")) |>
  rename(year = sampleYear,
         fish = Genetic.Sample.Number) |>
  filter(species == "Chinook") |>
  mutate(age = ageFresh + ageSalt +1) |>
  select(year, julian, fish, sexID, length, age)
  
  
#read in gsi and merge with asl-------------------------------------------------------------- 
gsi <- read.csv(here("analysis/data/raw/border-gsi-table.csv")) |>
  filter(year %in% c(2005:2019),
         gear == "Test Fishery") |>
  left_join( asl, by = c("year", "julian", "fish")) |>
  group_by(year, julian, fish) |>
  filter(prob == max(prob)) |>
  arrange(year, fish) |>
  as.data.frame()

age_comps <- gsi |>
  drop_na(age) |>
  group_by(year, CU, age) |>
  summarise(n = n()) |>
  arrange(as.numeric(age)) |>
  pivot_wider(names_from = age, 
              values_from = n) |>
  my_replace()|>
  arrange(year, CU) |>
  as.data.frame()


