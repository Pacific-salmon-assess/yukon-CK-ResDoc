## ---------- visualize HCRs -----------------------------------------------

out <- NULL

for(j in HCRs[c(2:6,18)]){
  HCR = j
  
  for(i in seq(1,500000, 500)) {
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
    if(grepl("fixed.ER.60", HCR)){
      ER <- ifelse(run.size==0, 0, 0.6)
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


# Visualize HCRs

ggplot(out) + geom_line(aes(x=run_size/1000, y=HR, col=HCR), linewidth=0.75) +
  facet_wrap(~factor(HCR, levels=HCRs[c(18,2:6)])) + theme_minimal() +
  scale_colour_manual(values=HCR_cols) +
  labs(x="Run Size (thousands)", y="Harvest Rate", title="Harvest Control Rules") +
  theme(legend.position = "none") +
  lims(x=c(0,400), y=c(0,1))
my.ggsave(here("analysis/plots/HCR_visualize.PNG"))

