### applying parameters ####
## rate and/or probabilities
#p = 10^-6
#x = p2/p1
#x = 1
#p2 = x*p/(x+1)
#p1 = p-p2
p1 = 10^-6
p2 = 10^-6
p3mr = 10^-5
p3cnd = 10^-5


# parameter for cells number for each genotype on doubling 0
# possible genotypes: "TKTK.TKTK"    "tk1TK.TKTK"   "tk2TK.TKTK"   "tk1tk1.TKTK"  "tk1TK.tk1TK"  "tk1tk2.TKTK" 
#"tk1TK.tk2TK"  "tk2TK.tk2TK"  "tk1tk1.tk1TK" "tk1tk1.tk2TK" "tk1tk2.tk1TK" "tk1tk2.tk2TK" "dead" "escaper"     
# two Cdels
cells <- c(1,rep(0,13))
# one Cdel
#cells <- c(0,0,0,1,rep(0,10))

# maximum number of doublings
max_cycle_cnt <- 70

### making transition matrix ####

# upload transition matrix
cdels2_transm <- read.csv("TKTKTKTK_transmatr_p3cnd.csv", stringsAsFactors=FALSE)
cdels2_transm_p3mr <- read.csv("TKTKTKTK_p3mr_transm.csv", stringsAsFactors=FALSE)

# calculating transition matrix p3mr
cdels2_transm_p3mr <- as.matrix(cdels2_transm_p3mr)
rownames(cdels2_transm_p3mr) <- cdels2_transm_p3mr[,1]
cdels2_transm_p3mr[,-1] -> cdels2_transm_p3mr

gsub("p3mr",p3mr,cdels2_transm_p3mr) -> cdels2_transm_p3mr_num

eval_all <- function(x) {
  return(as.numeric(eval(parse(text=x))))
}

apply(cdels2_transm_p3mr_num, c(1,2), eval_all) -> cdels2_transm_p3mr_num

# calculating transition matrix for others
cdels2_transm <- as.matrix(cdels2_transm)
rownames(cdels2_transm) <- cdels2_transm[,1]
cdels2_transm[,-1] -> cdels2_transm

gsub("p1",p1,cdels2_transm) -> cdels2_transm_num
gsub("p2",p2,cdels2_transm_num) -> cdels2_transm_num

eval_all <- function(x) {
  return(as.numeric(eval(parse(text=x))))
}

apply(cdels2_transm_num, c(1,2), eval_all) -> cdels2_transm_num

# check each row for reducibility to 1 
sum(cdels2_transm_num[1:14,])

## making all self-transitions 0 (not changing genotype cells are "leftowers" from the whole pool of cell particular genotype minus chenged one) 
for (i in 1:14) {
  cdels2_transm_num[i,i] <- 0
}


### writing function for Monte-Carlo simulations on the basis of this transmatrix ####

set.seed(42)

model_twoCdels <- function() {
  
  for(c in 1:max_cycle_cnt) {
    cells <- c(2*cells[1:12],cells[13],2*cells[14])
    
    t(t(cdels2_transm_num) * cells) -> lambda_p1p2p3cnd
    t(t(cdels2_transm_p3mr_num) * cells) -> lambda_p3mr
    
    tr_p1p2p3cnd <- apply(lambda_p1p2p3cnd, MARGIN = c(1,2), FUN = qpois, p=runif(1))
    tr_p3mr <- apply(lambda_p3mr, MARGIN = c(1,2), FUN = qpois, p=runif(1))
    
    tr_p3mr[1,2] -> tr_p3mr[4,2]
    tr_p3mr[1,3] -> tr_p3mr[13,3]
    tr_p3mr[2,5] -> tr_p3mr[9,5]
    tr_p3mr[4,6] -> tr_p3mr[13,6]
    tr_p3mr[2,7] -> tr_p3mr[13,7]
    tr_p3mr[3,7] -> tr_p3mr[10,7]
    tr_p3mr[3,8] -> tr_p3mr[13,8]
    tr_p3mr[4,9] -> tr_p3mr[14,9]
    tr_p3mr[4,10] -> tr_p3mr[13,10]
    tr_p3mr[6,11] -> tr_p3mr[14,11]
    tr_p3mr[9,11] -> tr_p3mr[13,11]
    tr_p3mr[6,12] -> tr_p3mr[13,12]
    tr_p3mr[10,12] + tr_p3mr[13,12] -> tr_p3mr[13,12]

    tr_p3mr + tr_p1p2p3cnd -> transitions
    
    cells - as.vector(colSums(transitions)) -> cells
    cells + as.vector(rowSums(transitions)) -> cells
    
     if (cells[length(cells)] > 0)
      return(c)
  }
  return(max_cycle_cnt + 1)
}


### calculations ####
# number of trials
trial_cnt <- 1000

# making Monte-Carlo simulations and writing them to the list with chosen /parameters/
MC_twoCdels <- replicate(trial_cnt, model_twoCdels())

# making a histogram out of it
hist(MC_twoCdels, freq = F)

# checking how many points we have below 1% quantile
sum(MC_twoCdels<quantile(MC_twoCdels, prob = 0.01))

# checking quantiles needed for FSL
quantile(MC_twoCdels, prob = c(0.01, 0.001, 0.0001, 0.00001))


### making doubling/FSL table out of MC simulations results ####
ceiling(quantile(MC_twoCdels, prob = c(0.10))) -> a
ceiling(quantile(MC_twoCdels, prob = c(0.05))) -> b
ceiling(quantile(MC_twoCdels, prob = c(0.01))) -> c
ceiling(quantile(MC_twoCdels, prob = c(0.001))) -> d
ceiling(quantile(MC_twoCdels, prob = c(0.0001))) -> e

if (sum(MC_twoCdels < a)<1 | sum(MC_twoCdels < b)<1 | sum(MC_twoCdels < c)<1 |
    sum(MC_twoCdels < d)<1 | sum(MC_twoCdels < e)<1) warning("Quantile is too low!!!")

as.data.frame(matrix(data=c(a, sum(MC_twoCdels < a), b, sum(MC_twoCdels < b),c, sum(MC_twoCdels < c),
                            d, sum(MC_twoCdels < d),e, sum(MC_twoCdels < e)), 
                     ncol=2, byrow = T)) -> MC_twoCdels_FSL
colnames(MC_twoCdels_FSL)<-c("doubling","success")
MC_twoCdels_FSL$trials <- rep(trial_cnt, nrow(MC_twoCdels_FSL))


# confidence intervals for MC
MC_twoCdels_FSL$upper <- rep(0,nrow(MC_twoCdels_FSL))
MC_twoCdels_FSL$lower <- rep(0,nrow(MC_twoCdels_FSL))


#install.packages("PropCIs")
library(PropCIs)

for (i in 1:nrow(MC_twoCdels_FSL)) {
  exactci(MC_twoCdels_FSL[i,2], MC_twoCdels_FSL[i,3], conf.level=0.95) -> a
  1/a$conf.int[1] -> MC_twoCdels_FSL[i,4]
  1/a$conf.int[2] -> MC_twoCdels_FSL[i,5]
}

MC_twoCdels_FSL$FSL <- MC_twoCdels_FSL$trials/MC_twoCdels_FSL$success


### nice plot with FSL ####
#install.packages("ggplot2")
library(ggplot2)
ggplot(MC_twoCdels_FSL, aes(x=log10(2^as.numeric(doubling)), y=log10(as.numeric(FSL))), les = 2) +
  geom_linerange(aes(ymax=log10(as.numeric(upper)), ymin=log10(as.numeric(lower))), size=0.5) +
  geom_point(size=1.5) +
  theme_light(base_size=15)+
  geom_line(size=0.5) +
  coord_cartesian(ylim=c(1,16), xlim=c(1,15)) +
  xlab("log10(cell number)") +
  ylab("log10(FSL)") +
  scale_y_continuous(breaks=c(seq(1:16))) +
  scale_x_continuous(breaks=c(seq(1:15)))


