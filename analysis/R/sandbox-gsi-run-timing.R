gsi <- read.csv(here("analysis/data/raw/border-gsi-table-2024-update.csv")) 

jds_2013 <- read.csv(here("analysis/data/raw/2013-gsi-jds.csv")) |>
  select(fish,JD,NewJD) |>
  distinct()
  
gsi_2013 <- gsi |>
  filter(year == 2013) |>
  left_join(jds_2013, by = "fish") |>
  mutate(julian = NewJD) |>
  select(-JD, -NewJD)

gsi_no_2013 <- gsi |>
  filter(year != 2013)

gsi_full <- rbind(gsi_no_2013,gsi_2013)

write.csv(gsi_full, here("analysis/data/raw/border-gsi-table-2024-update-full.csv"))


border_passage <- read.csv(here("analysis/data/raw/borderCounts.csv")) # border passage
gsi <- read.csv(here("analysis/data/raw/border-gsi-table-2024-update-full.csv")) 



# Figure S1: border passage and GSI samples -----
bc <- border_passage %>%
  drop_na(count) %>%
  group_by(year, count_type) %>%
  mutate(year_count = sum(count)) %>%
  mutate(max_count = max(count)) %>%
  ungroup() %>%
  mutate(prop = (count/year_count)) %>%
  mutate(prop2 = (count/max_count)) %>%
  mutate(year_gear = paste0(year,"_",count_type)) %>%
  filter(! year_gear %in% c("2005_fishWheel","2006_fishWheel","2007_fishWheel")) %>%
  filter(! year %in% c("1988","1989","1990","1998")) %>%
  as.data.frame()

gsi %>% 
  group_by(year, fish) %>%
  filter(prob == max(prob)) %>% 
  ungroup() %>%
  group_by(year) %>%
  summarise(med = median(prob)) %>% 
  ungroup() %>%
  mutate(max = max(med),
         min = min(med)) %>%
  summarise(quants = quantile(med, probs = c(0.05, 0.5, 0.95))) %>%
  as.data.frame()

gsi2 <- gsi %>%
  mutate(year_gear = paste0(year,"_",gear)) %>%
  group_by(year, fish) %>%
  filter(prob == max(prob)) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(year_count = n()) %>%
  ungroup() %>%
  group_by(year, julian) %>%
  mutate(julian_count = n()) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(max_count = max(julian_count)) %>%
  ungroup() %>%
  mutate(julian_prop = (julian_count/year_count)) %>%
  mutate(julian_prop2 = (julian_count/max_count)) %>%
  distinct(year, year_count, julian, julian_prop,julian_prop2) %>%
  as.data.frame()

gsi_count <- gsi |>
  group_by(year, fish) |>
  filter(prob == max(prob)) |>
  ungroup() |>
  group_by(year) |>
  summarize(year_count = n()) 
  
  
g <- ggplot() +
  geom_vline(xintercept=c(180,210,240), color="light grey",lwd=0.25,lty=2) +
  geom_bar(data = bc, aes(x=as.numeric(julian), y=prop2),stat = "identity") +
  geom_bar(data = gsi2, aes(x=as.numeric(julian), y=(julian_prop2)*-1), 
           fill="red", alpha=0.85, stat = "identity") +
  scale_x_continuous(limits=c(175,250), breaks = c(180,210,240)) +
  xlab("Day of year") +
  ylab("Run and GSI sample sizes") +
  facet_wrap(~year, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust = 1),
        strip.text = element_text(size=9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio=1,
        panel.spacing.x=unit(0.1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank()) +
  geom_text(data = gsi_count, 
            mapping = aes(x = 237, y = -0.5, label = year_count, hjust = 1, vjust = 2),
            size=3)

png(here("analysis/plots/trib-rr/gsi-run-samples.PNG"), width = 6, height = 8, units = "in", res = 600)
print(g)
dev.off()
