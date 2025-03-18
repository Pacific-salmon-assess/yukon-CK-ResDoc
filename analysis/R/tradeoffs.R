# function to find which element of x is closest to y
which_closest = function(x, y) {
  which.min(abs(x - y))
}

eq.ricker = function(ln_a, b) {
  
  Seq <- (ln_a-(-log(1-U.range)))/b
  Seq[Seq < 0] = 0

  Ceq <- ((ln_a-(-log(1-U.range)))/b)*exp(ln_a-b*((ln_a-(-log(1-U.range)))/b))-((ln_a-(-log(1-U.range)))/b)
  Ceq[is.na(Ceq)] = 0
  Ceq[Ceq < 0] = 0
  U_msy = (1 - lambert_W0(exp(1 - ln_a)))  
  overfished = ifelse(U.range > U_msy, 1, 0)
  extinct = ifelse(Seq == 0, 1, 0)
  
  return(list(S = Seq, C = Ceq, overfished = overfished, extinct = extinct))
  
}

par.posts.wide = par.posts |>
  pivot_wider(names_from = CU, values_from = c(ln_a, beta)) |>
  drop_na()

ni = length(par.posts.wide$sample)
ns = length(unique(post$CU))
U.range = seq(0, 1, 0.01)

n.U.range = length(U.range)

S.out = array(NA, dim = c(n.U.range,ni, ns))
C.out = array(NA, dim = c(n.U.range, ni, ns))
overfished.out = array(NA, dim = c(n.U.range, ni, ns))
extinct.out = array(NA, dim = c(n.U.range, ni, ns))

S.tot = matrix(NA, n.U.range, ni)
C.tot = matrix(NA, n.U.range, ni)
p.overfished = matrix(NA, n.U.range, ni)
p.extinct = matrix(NA, n.U.range, ni)

for (i in 1:ni) {
  for (s in 1:ns) {
    
    temp = eq.ricker(ln_a = as.numeric(par.posts.wide[i,s+1]) , b = as.numeric(par.posts.wide[i,s+10]) )
    
    S.out[,i,s] = temp$S
    C.out[,i,s] = temp$C
    overfished.out[,i,s] = temp$overfished
    extinct.out[,i,s] = temp$extinct
  }
  
  S.tot[,i] = rowSums(S.out[,i,])
  C.tot[,i] = rowSums(C.out[,i,])
  p.overfished[,i] = rowSums(overfished.out[,i,])/ns
  p.extinct[,i] = rowSums(extinct.out[,i,])/ns
}

S.tot.summ = apply(S.tot, 1, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.25, 0.75))))
C.tot.summ = apply(C.tot, 1, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.25, 0.75))))
overfished.summ = apply(p.overfished, 1, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.25, 0.75))))
extinct.summ = apply(p.extinct, 1, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.25, 0.75))))


tradeoffs <- cbind(U.range,S.tot.summ[3,],C.tot.summ[3,],overfished.summ[3,],extinct.summ[3,],
                   S.tot.summ[4,],C.tot.summ[4,],overfished.summ[4,],extinct.summ[4,],
                   S.tot.summ[5,],C.tot.summ[5,],overfished.summ[5,],extinct.summ[5,])
colnames(tradeoffs) <- c("U","spawn","harvest", "overfished", "extinct",
                         "spawn_lwr","harvest_lwr", "overfished_lwr", "extinct_lwr",
                         "spawn_upr","harvest_upr", "overfished_upr", "extinct_upr")
tradeoffs <- as.data.frame(tradeoffs)

# tradeoffs with aggregate exploitation rate 
ggplot(tradeoffs, aes(x=U*100, y=harvest/1000)) +
  geom_line(linetype = "solid", size=1.25, col = "grey") +
  geom_ribbon(aes(ymin=harvest_lwr/1000, ymax=harvest_upr/1000), alpha=0.2, lty=1,size=1.25, col = "grey") +
  geom_line(aes(x=U*100,y=overfished*50), lty=2, size=1.25, col = "orange") +
  geom_ribbon(aes(ymin=overfished_lwr*50, ymax=overfished_upr*50), lty=2,size=1.25, col = "orange", alpha=0.2, fill = "orange") +
  geom_line(aes(x=U*100,y=extinct*50), lty=3, size=1.25, col = "red") +
  geom_ribbon(aes(ymin=extinct_lwr*50, ymax=extinct_upr*50), lty=3,size=1.25, col = "red", alpha=0.2,fill = "red") +
  scale_y_continuous(sec.axis = sec_axis(~. * 2, name = "Populations at risk (%)")) +
  ylab("Harvest (000s)") +
  xlab("Harvest rate (%)") +
  theme_bw() +
  theme(axis.title = element_text(size= 9),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,0,3,0.5), units = "lines")) +
  annotate("text", x = c(19,19,19),
           y = c(45,40,35),
           label = c("Harvest", "Overfished", "Extinct"),
           color=c("grey","orange","red"), 
           size=3.5,
           hjust=0) +
  annotate("segment", x = c(14,14,14),
           xend=c(17,17,17),
           y = c(45,40,35),
           yend = c(45,40,35),
           lty = c(1,2,3),
           color=c("grey","orange","red"), 
           size=1)
my.ggsave(here("analysis/plots/U-tradeoffs.PNG"))


# tradeoffs with aggregate escapement 
ggplot(tradeoffs, aes(x=spawn/1000, y=harvest/1000)) +
  geom_line(linetype = "solid", size=1.25, col = "grey") +
  geom_line(aes(x=spawn/1000,y=overfished*50), lty=2, size=1.25, col = "orange") +
  geom_line(aes(x=spawn/1000,y=extinct*50), lty=3, size=1.25, col = "red") +
  scale_y_continuous(sec.axis = sec_axis(~. * 2, name = "Populations at risk (%)")) +
  ylab("Harvest (000s)") +
  xlab("Spawners (000s)") +
  theme_bw() +
  theme(axis.title = element_text(size= 9),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,0,3,0.5), units = "lines")) +
  annotate("text", x = c(79,79,79),
           y = c(45,40,35),
           label = c("Harvest", "Overfished", "Extinct"),
           color=c("grey","orange","red"), 
           size=3.5,
           hjust=0) +
  annotate("segment", x = c(74,74,74),
           xend=c(77,77,77),
           y = c(45,40,35),
           yend = c(45,40,35),
           lty = c(1,2,3),
           color=c("grey","orange","red"), 
           size=1)


ggplot(tradeoffs, aes(x=spawn/1000, y=harvest/1000)) +
  geom_line(linetype = "solid", size=1.25, col = "grey") +
  geom_ribbon(aes(ymin=harvest_lwr/1000, ymax=harvest_upr/1000), alpha=0.2, lty=1,size=1.25, col = "grey") +
  geom_line(aes(x=spawn/1000,y=overfished*50), lty=2, size=1.25, col = "orange") +
  geom_ribbon(aes(ymin=overfished_lwr*50, ymax=overfished_upr*50), lty=2,size=1.25, col = "orange", alpha=0.2, fill = "orange") +
  geom_line(aes(x=spawn/1000,y=extinct*50), lty=3, size=1.25, col = "red") +
  geom_ribbon(aes(ymin=extinct_lwr*50, ymax=extinct_upr*50), lty=3,size=1.25, col = "red", alpha=0.2,fill = "red") +
  scale_y_continuous(sec.axis = sec_axis(~. * 2, name = "Populations at risk (%)")) +
  ylab("Harvest (000s)") +
  xlab("Spawners (000s)") +
  theme_bw() +
  theme(axis.title = element_text(size= 9),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,0,3,0.5), units = "lines")) +
  annotate("text", x = c(79,79,79),
           y = c(45,40,35),
           label = c("Harvest", "Overfished", "Extinct"),
           color=c("grey","orange","red"), 
           size=3.5,
           hjust=0) +
  annotate("segment", x = c(74,74,74),
           xend=c(77,77,77),
           y = c(45,40,35),
           yend = c(45,40,35),
           lty = c(1,2,3),
           color=c("grey","orange","red"), 
           size=1)

my.ggsave(here("analysis/plots/s-tradeoffs.PNG"))
