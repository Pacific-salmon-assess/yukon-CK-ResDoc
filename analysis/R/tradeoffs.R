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

post = bench.par.table |>
  filter(bench.par %in% c("alpha","beta")) |>
  select(CU, bench.par, mean) |>
  pivot_wider(names_from = bench.par, values_from = mean) |>
  mutate(ln_alpha = log(alpha))

ni = 1
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
    
    temp = eq.ricker(ln_a = as.numeric(post[s,4]), b = as.numeric(post[s,3]))
    
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

plot(U.range,S.tot, type="l")
lines(U.range,C.tot, type="l")

S.tot.summ = apply(S.tot, 1, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.025, 0.975))))
C.tot.summ = apply(C.tot, 1, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.025, 0.975))))
overfished.summ = apply(p.overfished, 1, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.025, 0.975))))
extinct.summ = apply(p.extinct, 1, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.025, 0.975))))


tradeoffs <- cbind(U.range,S.tot,C.tot,overfished.summ[1,],extinct.summ[1,])
colnames(tradeoffs) <- c("U","spawn","harvest", "overfished", "extinct")
tradeoffs <- as.data.frame(tradeoffs)

ggplot(tradeoffs, aes(x=U*100, y=harvest/500)) +
  geom_line(linetype = "solid", size=0.5) +
  geom_line(aes(x=U*100,y=overfished*100), lty=2, size=0.5) +
  geom_line(aes(x=U*100,y=extinct*100), lty=3, size=0.5) +
  scale_y_continuous(limits = c(0,30), breaks=c(0,10,20,30,40)) +
  scale_y_continuous(sec.axis = sec_axis(~./1, name = "Populations at risk (%)"),limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100), breaks=c(0,20,40,60,80,100)) +
  ylab("Harvest (000s)") +
  xlab("Harvest rate (%)") +
  theme_bw() +
  theme(axis.title = element_text(size= 9),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,0,3,0.5), units = "lines")) +
  annotate("text", x = c(9,9,9),
           y = c(95,90,85),
           label = c("Harvest", "Overfished", "Extinct"),
           color="black", 
           size=2,
           hjust=0) +
  annotate("segment", x = c(0,0,0),
           xend=c(7,7,7),
           y = c(95,90,85),
           yend = c(95,90,85),
           lty = c(1,2,3),
           color="black", 
           size=0.3)
