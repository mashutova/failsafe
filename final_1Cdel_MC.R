### applying parameters ####
## probabilities and rate
#p = 10^-6
# x = p2/p1
#x = 1
#p2 = x*p/(x+1)
#p1 = p-p2
p1=10^-6
p2=10^-6
p3mr=10^-5
p3cnd=10^-5

# maximum number of doublings
max_cycle_cnt <- 50

# parameter for cells number for each genotype on doubling 0
# possible genotypes: "TKTK" "TKtk1" "TKtk2" "dead" "escaper"    
cells <- c(1,rep(0,4))
#cells <- c(0,1,0,0,0)

### making transition matrix ####
# upload transition matrix
cdel1_transm <- read.csv("/Users/mashut/Documents/Nagy_lab/failsafe/finalised_failsafe_scripts/failsafe_clean/final_transmatrixes/TKTK_transmatr_p1p2p3cndp3mr_1.csv", stringsAsFactors=FALSE)

# calculating transition matrix
cdel1_transm <- as.matrix(cdel1_transm)
rownames(cdel1_transm) <- cdel1_transm[,1]
cdel1_transm[,-1] -> cdel1_transm

cdel1_transm[1,2:3] <- c(0,0)

gsub("p1",p1,cdel1_transm) -> cdel1_transm_num
gsub("p2",p2,cdel1_transm_num) -> cdel1_transm_num

eval_all <- function(x) {
  return(as.numeric(eval(parse(text=x))))
}

apply(cdel1_transm_num, c(1,2), eval_all) -> cdel1_transm_num

# check each row for reducibility to 1 
sum(cdel1_transm_num[1:5,])

## making all self-transitions 0 (not changing genotype cells are "leftowers" from the whole pool of cell particular genotype minus chenged one) 
for (i in 1:5) {
  cdel1_transm_num[i,i] <- 0
}


### writing function for Monte-Carlo simulations on the basis of this transmatrix ####
set.seed(42)

model_oneCdel_1 <- function() {
  
  for(c in 1:max_cycle_cnt) {
    cells <- c(2*cells[1:3],cells[4],2*cells[5])
    
    t(t(cdel1_transm_num) * cells) -> lambda_p1p2p3
    transitions <- apply(lambda_p1p2p3, MARGIN = c(1,2), FUN = qpois, p=runif(1))
    
    cells - as.vector(colSums(transitions)) -> cells
    cells + as.vector(rowSums(transitions)) -> cells
    
    if (cells[length(cells)] > 0)
      return(c)
  }
  return(max_cycle_cnt + 1)
}


### calculations ####
# number of trials
trial_cnt <- 10000

# making Monte-Carlo simulations and writing them to the list with chosen /parameters/
MC_oneCdel1 <- replicate(trial_cnt, model_oneCdel_1())

# making a histogram out of it
hist(MC_oneCdel, freq = F)

# checking how many points we have below 1% quantile
sum(MC_oneCdel<quantile(MC_oneCdel, prob = 0.01))

# checking quantiles needed for FSL
quantile(MC_oneCdel, prob = c(0.1, 0.01, 0.001, 0.0001, 0.00001))


### making doubling/FSL table out of MC simulations results ####
ceiling(quantile(MC_oneCdel, prob = c(0.10))) -> a
ceiling(quantile(MC_oneCdel, prob = c(0.01))) -> b
ceiling(quantile(MC_oneCdel, prob = c(0.001))) -> c
ceiling(quantile(MC_oneCdel, prob = c(0.0001))) -> d
ceiling(quantile(MC_oneCdel, prob = c(0.00001))) -> e

if (sum(MC_oneCdel < a)<1 | sum(MC_oneCdel < b)<1 | sum(MC_oneCdel < c)<1 |
    sum(MC_oneCdel < d)<1 | sum(MC_oneCdel < e)<1) warning("quantile is too low!!!")

as.data.frame(matrix(data=c(a, sum(MC_oneCdel < a), b, sum(MC_oneCdel < b),c, sum(MC_oneCdel < c),
                            d, sum(MC_oneCdel < d),e, sum(MC_oneCdel < e)), 
                     ncol=2, byrow = T)) -> MC_oneCdel_FSL
colnames(MC_oneCdel_FSL)<-c("doubling","success")
MC_oneCdel_FSL$trials <- rep(trial_cnt, nrow(MC_oneCdel_FSL))


# confidence intervals for MC
MC_oneCdel_FSL$upper <- rep(0,nrow(MC_oneCdel_FSL))
MC_oneCdel_FSL$lower <- rep(0,nrow(MC_oneCdel_FSL))


#install.packages("PropCIs")
library(PropCIs)

for (i in 1:nrow(MC_oneCdel_FSL)) {
  exactci(MC_oneCdel_FSL[i,2], MC_oneCdel_FSL[i,3], conf.level=0.95) -> a
  1/a$conf.int[1] -> MC_oneCdel_FSL[i,4]
  1/a$conf.int[2] -> MC_oneCdel_FSL[i,5]
}

MC_oneCdel_FSL$FSL <- MC_oneCdel_FSL$trials/MC_oneCdel_FSL$success


### nice plot with FSL ####
#install.packages("ggplot2")
library(ggplot2)
ggplot(MC_oneCdel_FSL, aes(x=log10(2^as.numeric(doubling)), y=log10(as.numeric(FSL))), les = 2) +
  geom_linerange(aes(ymax=log10(as.numeric(upper)), ymin=log10(as.numeric(lower))), size=0.5) +
  #geom_boxplot(aes(fill=treatment), alpha = 0.2) +
  geom_point(size=1.5) +
  #facet_wrap(~ dose) +
  #scale_colour_hue(l=50) + # Use a slightly darker palette than normal
  #geom_smooth(se=F, color="black", aes(group=1)) +
  theme_light(base_size=15)+
  geom_line(size=0.5) +
  coord_cartesian(ylim=c(1,16), xlim=c(1,15)) +
  xlab("log10(cell number)") +
  ylab("log10(FSL)") +
  scale_y_continuous(breaks=c(seq(1:16))) +
  scale_x_continuous(breaks=c(seq(1:15)))

