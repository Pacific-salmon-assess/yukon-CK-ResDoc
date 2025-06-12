#wrangling ASL data i

library(tidyverse)
library(here)
library(zoo)
library(ggpubr)
library(ggsidekick)

# Hamachan selectivity corrected age and sex comps 1981-2006
fw_age_sex <- read.csv(here("analysis/data/raw/fw-border-age-sex-comps.csv"))

# Hamachan length raw data with selectivity 1981-2006
fw_l <- read.csv(here("analysis/data/raw/FWborder.csv"))

# Eage ASL data 2005-2024
eagle_age_sex_len <- read.csv(here("analysis/data/raw/ASL_Output_Chinook_Eagle_2005-2024.csv"))

## wrangling for full age and sex dataset
#  eagle ASL prep
eagle_asl <- eagle_age_sex_len %>% 
  filter(!is.na(Total.Age), Sex != "unknown", Species != "Chum") |>
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
write.csv(cdn_as_1985_2024[which(cdn_as_1985_2024$Sex=="male"),], here("analysis/data/raw/male_age_comps.csv"),row.names = FALSE)

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
  filter(!is.na(Total.Age), Sex != "unknown", Species != "Chum") |>
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


#### plots

fem_len_comp <- read.csv(here("analysis/data/raw/female_length_comps.csv"))
fem_age_comps <- read.csv(here("analysis/data/raw/female_age_comps.csv"))
fem_egg_mass_comps <- read.csv(here("analysis/data/raw/female_length_comps_egg_mass.csv")) # female spawner eggs by age
fem_egg_comps <- read.csv(here("analysis/data/raw/female_length_comps_eggs.csv")) # female spawner eggs by age

reproOutputPerSpawner <- rowSums(fem_S_comps*fem_egg_mass_comps)
reproOutputPerSpawnerEggs <- rowSums(fem_S_comps*fem_egg_comps)
reproOutput <- data.frame(cbind(seq(1985,2024),reproOutputPerSpawnerEggs))
reproOutputEM <- data.frame(cbind(seq(1985,2024),reproOutputPerSpawner))
write.csv(reproOutputEM, here("analysis/data/generated/reproOutputEM.csv"),row.names = FALSE)

fem_age_comp <- fem_age_comps[,c(1,3:6)]%>%
  group_by(Sample.Year) %>%
  mutate(total_prop=sum(age_4,age_5,age_6,age_7),
         age4=age_4/total_prop,
         age5=age_5/total_prop,
         age6=age_6/total_prop,
         age7=age_7/total_prop)%>%
  select(Sample.Year,age4, age5, age6, age7)%>%
  rename(Four = age4, Five = age5, Six = age6, Seven = age7)%>%
  pivot_longer(!Sample.Year , names_to = "Age",
               values_to = "prop")

fem_age_comp$age_f <- factor(fem_age_comp$Age, levels = c("Four", "Five", "Six", "Seven"))

a <- ggplot(fem_age_comp, aes(fill=age_f, y=prop, x=Sample.Year)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Year") +
  ylab("Proportion") +
  scale_fill_viridis_d(name = "Age") +
  theme_sleek() +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.position="top",
        plot.margin = margin(10,20,0.5,20))


prop_females <- rowSums(fem_age_comps[,3:6])

sex_ratio<-cbind(seq(1985,2024),prop_females)
colnames(sex_ratio)<-c("Year","prop_fem")
sex_ratio<-as.data.frame(sex_ratio)

b <- ggplot(sex_ratio, aes(x = Year, y = prop_fem)) +
  geom_smooth(method="lm", color="grey") +
  geom_point(size=2, color="dark grey")+ 
  xlab("Year") +
  ylab("Proportion female") +
  coord_cartesian(ylim=c(0,1)) +
  theme_sleek() +
  theme(plot.margin = margin(15,20,0.5,20))

laa <- as.data.frame(fem_len_comp)%>%
  select(Sample.Year,age_4, age_5, age_6, age_7)%>%
  rename(Four = age_4, Five = age_5, Six = age_6, Seven = age_7)%>%
  pivot_longer(!Sample.Year , names_to = "Age",
               values_to = "length") 


laa$age_f <- factor(laa$Age, levels = c("Four", "Five", "Six", "Seven"))

c <- ggplot(laa, aes(x = Sample.Year, y = length)) +
  geom_smooth(method="lm", color="grey") +
  geom_point(size=2, color="dark grey")+ 
  xlab("Year") +
  ylab("Female length \n  (mm; MEFL)") +
  theme_sleek() +
  facet_wrap(~age_f, scales = "free_y") +
  scale_x_continuous(breaks=c(1985, 1995, 2005, 2015)) +
  theme(strip.text.x = element_text(size=8),
        plot.margin = margin(0.5,10,1,0.5))

d <- ggplot(reproOutput, aes(x = V1, y = reproOutputPerSpawnerEggs)) +
  geom_smooth(method="lm", color="grey") +
  geom_point(size=2, color="dark grey")+ 
  xlab("Year") +
  ylab("Average reproductive output \n (total eggs per spawner)") +
  theme_sleek() +
  theme(plot.margin = margin(45,10,0.5,0.5))
g <- ggarrange(b,c,a,d, 
               labels = c("a", "b","c", "d"),
               heights = c(0.8,1))


g
my.ggsave(here("analysis/plots/trib-rr/asl.PNG"))
ggsave(here("csasdown/figure/asl.PNG"), height = 550*2, 
       width = 700*2, units="px", dpi=200)


# calculate total egg and egg mass time series for Jim M
rr <- read.csv(here("analysis/data/raw/rr-table.csv"))

cdn_rr <- rr %>%
  filter(stock == "Canadian",
         Year >1984) |>
  select(Year,Escapement)

cdn_rr$total_females <- cdn_rr$Escapement*rowSums(fem_S_comps)
cdn_rr$total_eggs <- rowSums(cdn_rr$Escapement*fem_S_comps*female_length_comps_eggs)
cdn_rr$total_egg_mass <- rowSums(cdn_rr$Escapement*fem_S_comps*female_length_comps_egg_mass)

write.csv(cdn_rr, here("analysis/data/raw/cdn-yukon-chinook-repro-output.17Mar2025.csv"),row.names = FALSE)

# age by CU ----

eagle_age_sex_gen <- read.csv(here("analysis/data/raw/ASL_Eagle_2005-2024_geneticIDs.csv"))
gsi <- read.csv(here("analysis/data/raw/border-gsi-table-2024-update-full.csv")) 

age <-eagle_age_sex_gen |>
  mutate(fish=Genetic.Sample.Number,
         year=sampleYear) |>
  filter(species == "Chinook",
         year>2008) |>
  select(year, fish, sexID, totalAge)

gsi_hProb <- gsi |>
  filter(prob>0.5,
         year>2008) |>
  select(year, fish, CU, prob)

gsi_age <- age |>
  left_join(gsi_hProb, by = c("year", "fish")) |>
  mutate(age = case_when(totalAge == 1.1 ~ 3,
                         totalAge == 1.2 ~ 4,
                         totalAge == 1.3 ~ 5,
                         totalAge == 1.4 ~ 6,
                         totalAge == 1.5 ~ 7,
                         totalAge == 2.2 ~ 5,
                         totalAge == 2.3 ~ 6,
                         totalAge == 2.4 ~ 7)) |>
  filter(age != 3)


cu_age_sex <- gsi_age |>
  group_by(year,CU,sexID) |>
  count(age) |>
  drop_na()

cu_age <- gsi_age |>
  group_by(year,CU) |>
  count(age) |>
  drop_na()|>
  mutate(year_count = sum(n),
         prop = n/year_count) |>
  drop_na()

cu_age_all_yr <- gsi_age |>
  group_by(CU) |>
  count(age) |>
  mutate(year_count = sum(n),
         prop = n/year_count) |>
  drop_na()

ggplot(cu_age_all_yr, aes(x = age, y = prop)) +
  geom_bar(position="dodge", stat = "identity") + 
  theme_sleek() +
  facet_wrap(~CU) 

ggplot(cu_age |> filter(!year %in% c(2010,2012,2013)), aes(x = CU, y = prop, fill=as.factor(age))) +
  geom_bar( stat = "identity") + 
  theme_sleek() +
  facet_wrap(~year) +
  theme(legend.position = c(0.7,0.145),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Conservation Unit", y = "Proportion") 
my.ggsave(here("analysis/plots/SR_models/age-cu-by-yrs.PNG"))

ggplot(cu_age_all_yr, aes(x = CU, y = prop, fill=as.factor(age))) +
  geom_bar( stat = "identity") + 
  theme_sleek() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Conservation Unit", y = "Proportion") 
my.ggsave(here("analysis/plots/SR_models/age-cu-all-yrs.PNG"))
