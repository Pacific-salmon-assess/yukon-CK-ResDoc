## ---------- visualize HCRs -----------------------------------------------

out <- NULL

for(j in HCRs[2:6]){
  HCR = j
  
  for(i in seq(1,500000, 1000)) {
    run.size = i
    
    if(is.na(run.size)==TRUE){run.size <- 0}
    if(run.size > 999000) {run.size <- 1000000} 
    if(HCR == "no.fishing"){HR.all <- 0}
    if(HCR == "status.quo"){ 
      catch <- ifelse(run.size<=42500, 0, run.size-42500)
      HR.all <- ifelse(run.size==0, 0, catch/run.size)
      if(HR.all > 0.8){       ## Add ER cap (80%) 
        catch <- run.size*0.8
        HR.all <- catch/run.size }}
    if(HCR == "status.quo.cap"){
      catch <- ifelse(run.size<=42500, 0, run.size-42500)
      HR.all <- ifelse(run.size==0, 0, catch/run.size)
      if(HR.all > 0.4){       ## Lower ER cap (40%) 
        catch <- run.size*0.4
        HR.all <- catch/run.size }}
    if(HCR == "rebuilding"){
      catch <- ifelse(run.size<=71000, 0, run.size-71000)
      HR.all <- ifelse(run.size==0, 0, catch/run.size)
      if(HR.all > 0.8){       ## ER cap (80%) 
        catch <- run.size*0.8
        HR.all <- catch/run.size }}
    if(HCR == "rebuilding.cap"){
      catch <- ifelse(run.size<=71000, 0, run.size-71000)
      HR.all <- ifelse(run.size==0, 0, catch/run.size)
      if(HR.all > 0.4){       ## lower ER cap (40%) 
        catch <- run.size*0.4
        HR.all <- catch/run.size }}
    if(grepl("fixed.ER", HCR)){
      ER <- ifelse(run.size==0, 0, as.numeric(gsub("\\D", "", HCR)))
      catch <- run.size*ER
      HR.all <- ER}
    if(HCR == "alt.rebuilding"){
      if(run.size <= 19000) catch <- 0
      if(run.size >= 95000) catch <- run.size*0.4
      if(run.size > 19000 & run.size < 95000){
        dat <- data.frame(R=c(19000,95000), ER=c(0,0.4))
        lin <- lm(ER ~ R, data=dat)
        ER <- coef(lin)[1] + coef(lin)[2]*run.size
        catch <- run.size*ER
      }
      HR.all <- catch/run.size
    }
    itn <- data.frame(HCR = HCR,
                      run_size=run.size,
                      HR = HR.all,
                      catch=catch)
    out <- bind_rows(out, itn)
  }
}

# Add fixed ER to out
out_fixed_er <- data.frame(HCR="fixed.ER",
                           run_size=rep(seq(1,500000, 1000), times=20),
                           HR = rep(seq(0.05,1,0.1), each=20),
                           ER = rep(seq(0.05,1,0.1), each=20),
                           catch = NA) 
out_plot <- bind_rows(out, out_fixed_er)

# load historical run size info
hist_run <- read.csv(here('analysis', 'data', 'raw', 'rr-table.csv'))
hist_run <- hist_run %>% dplyr::summarize(lower = quantile(Total.run, 0.025), 
                                          upper = quantile(Total.run, 0.975))

# Visualize HCRs

ggplot(out) + geom_line(aes(x=run_size/1000, y=HR*100, col=HCR), linewidth=0.75) +
  geom_rect(data=hist_run, aes(xmin = lower/1000, xmax=upper/1000, ymin=0, ymax=100), fill="grey70", alpha=0.2) +
  facet_wrap(~factor(HCR, levels=unique(out$HCR)[c(3:5,1:2)])) + 
  scale_colour_manual(values=HCR_cols) +
  labs(x="Run Size (thousands)", y="Harvest Rate (%)") +
  theme_minimal() + theme(legend.position = "none") 
  lims(x=c(0,400)) +
  scale_y_continuous(breaks=seq(0,100,20), limits=c(0,100)) 
my.ggsave(here("analysis/plots/HCR_visualize.PNG"))

