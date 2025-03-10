bench.par1 <- read.csv(here("analysis/data/generated/bench_par_table.csv"))
bench.par2 <- read.csv(here("analysis/data/generated/brendan_bench_par_table.csv"))
bench.par3 <- read.csv(here("analysis/data/generated/hannah_bench_par_table.csv"))

#Want all these to be TRUE, run inside of all() to see WHICH dont match
all(bench.par1$CU == bench.par2$CU)
#if above is FALSE, which CUs are getting mixed up?
cbind(bench.par1$CU, bench.par2$CU)[which(bench.par1$CU != bench.par2$CU),]

all(bench.par1$bench.par == bench.par2$bench.par)
all(bench.par1$mean == bench.par2$mean) #betas match

cbind(bench.par1$CU, bench.par2$CU)
cbind(bench.par1$CU, bench.par3$CU)

# join benchmark tables by CU 
bench.join <- left_join(bench.par1[,1:3], bench.par2[,1:3], by=c("CU", "bench.par"))
bench.join <- bench.join %>% mutate(diff=abs(mean.x - mean.y)) %>% mutate(diff_mag = (diff/mean.x)*100)
View(bench.join %>% group_by(bench.par) %>% arrange(desc(diff), .by_group = T))

# visualize
for(bench in unique(bench.join$bench.par)[1:7]){
  bench.join.sub <- filter(bench.join, bench.par %in% unique(bench.join$bench.par)[1:6])
  
    ggplot(bench.join.sub) + geom_point(aes(x=CU, y=diff_mag, col=bench.par), size=2) + 
    geom_hline(aes(yintercept=mean(diff_mag))) +
    labs(title="benchmark error", x="", y="benchmark diff / benchmark value") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle=0))
  

  pars.join.sub <- filter(bench.join, bench.par %in% unique(bench.join$bench.par)[7:10])
  
    ggplot(pars.join.sub) + geom_point(aes(x=CU, y=diff_mag, col=bench.par), size=2) + 
    geom_hline(aes(yintercept=mean(diff_mag))) +
    labs(title = "parameter error", x="", y="param diff / param value") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle=0))
  
}


#look at difference in means
diffs <- data.frame(CU = bench.par1$CU, 
                    bench.par = bench.par1$bench.par,
                    diff = as.numeric(bench.par1$mean - bench.par2$mean))
diffs %>% summarize(highest=max(diff), .by="bench.par")
View(diffs %>% group_by(bench.par) %>% arrange(desc(diff), .by_group = T))

cbind(bench.par1$CU, bench.par2$CU)

#is it just the order that is bessed up? sort means and compare
cbind(sort(bench.par1$mean), sort(bench.par2$mean))
