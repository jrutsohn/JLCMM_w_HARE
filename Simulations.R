#Note: Simulations shown here reflect original submission.  New simulations have been added (Replicated_Simulations), and some simulations presented
## below have been placed in the online companion (09-01-2024).
library(tidyverse)
library(dplyr)
library(ggplot2)
library(magrittr)
library(lcmm)
library(polspline)
library(survival)
library(mvtnorm)

############################################
#Producing the design matrix for covariates
#
#Using named variables for readability
# Age is x1 in manuscript
# Trt is x2 in manuscript
# Sex is x3 in manuscript
############################################
design_matrix <- function(seed,n){
  set.seed(seed)
  #Sex: Men (sex = 0), Women (sex = 1)
  sex <- ifelse(runif(n) < 0.5, 0, 1)
  #Age (scaled and centered)
  age <- c(runif(n*0.2, min = 30, max = 65), runif(n*0.4, min = 65, max = 75), runif(n*0.40, min = 75, max = 85))
  age <- (age - 70)/10
  trt <- ifelse(runif(n) < 0.5, 0, 1)
  X <- cbind(sex,age,trt)
  colnames(X) <- c("sex","age","trt")
  return(X)
}

#Simulating design matrices for three latent classes
X1 <- design_matrix(2702,200)
X2 <- design_matrix(2702,150)
X3 <- design_matrix(2702,50)

#Regression parameters for the three latent classes.  See Table 1.
beta1 <- c(1.5, 0.5, 0.25)
beta2 <- c(0.5, -0.5, -0.25)
beta3 <- c(3.0, 1.0, 0.125)


gamma1 <- c(1.125, 0.25)
gamma2 <- c(-0.05, 0.55)
gamma3 <- c(0.35, -0.35)


#Scale
lambda1 <- 1
lambda2 <- 0.9
lambda3 <- 1.1

#Shape
psi1 <- 0.5
psi2 <- 1.5
psi3 <- 3


#################################
#Function that simulates
#longitudinal data and
#survival data based on
#design matrix & parameters
#################################

design_all <- function(X,beta,lambda,psi,gamma,IDoffset){
  #X: Design matrix of covariates
  #beta: coefficients for mixed-effects model
  #IDoffset: Number to start ID (to make unique cases)
  n <- dim(X)[1]
  ID <- seq(1,n,1) + IDoffset
  times <- sample(1:5, n, replace = T)
  visits <- c()
  for(i in 1:n){
    visits[[i]] <- seq(1, times[i], 1)/2 - 0.5 #Sets minimum to 0 and scales from 0 to 2
  }
  vv <- data.frame()
  for(j in 1:n){
    a <- data.frame((cbind(ID[j],visits[[j]])))
    vv <- rbind(vv, a)
  }
  
  #sigma.2: Matrix of random effects where v1 = 0.4, v2 = 0.2, rho = 0.4
  sigma.2 <- rbind(c(0.4^2, 0.2*0.4^2), c(0.2*0.4^2, 0.2^2))
  bi <- rmvnorm(n, mean = rep(0, 2), sigma = sigma.2)
  e <- rnorm(dim(vv)[1], 0, 0.1) #Note: changed epsilon for scaling reasons
  
  des <- cbind(ID,X,data.frame(bi))
  names(des) <- c("ID","sex","age","trt","b0","b1")
  err <- cbind(vv,e)
  names(err) <- c("ID","visit","epsilon")
  
  mod <- left_join(err,des,by="ID") %>% select(ID,visit,sex,age,trt,b0,b1,epsilon)
  
  surv.inv.mu <- mod %>% mutate(mu = beta[1] + beta[2]*visit + beta[3]*age + 0.35*trt + 0.5*sex + b0 + b1*visit) %>%
    group_by(ID) %>% filter(row_number() == n()) %>% ungroup()
  longout <- mod %>% mutate(y = beta[1] + beta[2]*visit + beta[3]*age + 0.35*trt + 0.5*sex + b0 + b1*visit + epsilon)
  
  surv.mu <- surv.inv.mu$mu
  
  
  u <- runif(n,0,1)
  xtb <- gamma[1]*des$trt + gamma[2]*des$age + 0.05*surv.mu #0.05 association coefficient
  s.inv <- (-log(u)/(lambda*exp(xtb)))^(1/psi)
  delta <- as.vector(s.inv < 2.1)
  s.inv <- as.vector(ifelse(delta,s.inv,2.1))
  delta <- as.numeric(delta)
  
  surv.times <- data.frame(cbind(ID,s.inv,delta))
  dat <- left_join(longout,surv.times,by="ID")
  
  return(dat)
}

Y1 <- design_all(X1, beta = beta1, lambda = lambda1, psi = psi1, gamma = gamma1, IDoffset = 0) #IDs 0 - 200
Y2 <- design_all(X2, beta = beta2, lambda = lambda2, psi = psi2, gamma = gamma2, IDoffset = 200) #IDs 201 - 350
Y3 <- design_all(X3, beta = beta3, lambda = lambda3, psi = psi3, gamma = gamma3, IDoffset = 350) #IDs 351 - 400

#Combinint the datasets
YS <- data.frame(rbind(Y1,Y2,Y3))

#Covariate list for HARE
cov.list <- c("sex","age","trt")


#HARE model
hare.mod <- hare(data = YS$s.inv, delta = YS$delta, cov = YS[,cov.list]) #Specifies the time to event (s.inv), event indicator (delta = {0,1}), and covariate list (YS[,cov.list])

#Extracting knots from HARE model
hare.mod$knots

#Knot numbers and sequence to insert into JLCMMs for HARE-based models
knot.k <- hare.mod$knots[1,1]
hazard.k <- knot.k + 2
hazard.str <- paste0(hazard.k,"-manual-splines")
knot.seq <- c(hare.mod$knots[1,-1])


###########################################
#JLCMMs using Jointlcmm() function
# from lcmm package.
#
#Note that gridsearch() is used for HARE
# to ensure finding solution that reaches
# global maximum.
##########################################

weibull1 <- Jointlcmm(y ~ visit + age + sex + trt,
                      random =~ visit,
                      survival = Surv(s.inv, delta) ~ age + sex + trt,
                      hazard = "Weibull",
                      logscale = FALSE,
                      subject = 'ID',
                      data = YS,
                      ng = 1)

hare1 <- Jointlcmm(y ~ visit + age + sex + trt,
                   random =~ visit,
                   survival = Surv(s.inv, delta) ~ age + trt,
                   hazard = hazard.str,
                   hazardnodes = knot.seq[1:knot.k],
                   logscale = FALSE,
                   subject = 'ID',
                   data = YS,
                   ng = 1)

weibull2 <- Jointlcmm(y ~ visit + age + sex + trt,
                      random =~ visit,
                      mixture =~ age + visit,
                      survival = Surv(s.inv, delta) ~ mixture(age) + mixture(sex) + mixture(trt),
                      hazard = "Weibull",
                      logscale = FALSE,
                      subject = 'ID',
                      data = YS,
                      B = weibull1,
                      ng = 2)

weibull3 <- Jointlcmm(y ~ visit + age + sex + trt,
                      random =~ visit,
                      mixture =~ age + visit,
                      survival = Surv(s.inv, delta) ~ mixture(age) + mixture(sex) + mixture(trt),
                      hazard = "Weibull",
                      logscale = FALSE,
                      subject = 'ID',
                      data = YS,
                      B = weibull1,
                      ng = 3)

weibull4 <- Jointlcmm(y ~ visit + age + sex + trt,
                      random =~ visit,
                      mixture =~ age + visit,
                      survival = Surv(s.inv, delta) ~ mixture(age) + mixture(sex) + mixture(trt),
                      hazard = "Weibull",
                      logscale = FALSE,
                      subject = 'ID',
                      data = YS,
                      B = weibull1,
                      ng = 4)

hare2 <- gridsearch(rep = 30, maxiter = 15, minit = hare1,
                    Jointlcmm(y ~ visit + age + sex + trt,
                              random =~ visit,
                              mixture =~ age + visit,
                              survival = Surv(s.inv, delta) ~ mixture(age) + mixture(trt),
                              hazard = hazard.str,
                              hazardnodes = knot.seq[1:knot.k],
                              logscale = FALSE,
                              subject = 'ID',
                              data = YS,
                              ng = 2))

hare3 <- gridsearch(rep = 30, maxiter = 15, minit = hare1,
                    Jointlcmm(y ~ visit + age + sex + trt,
                              random =~ visit,
                              mixture =~ age + visit,
                              survival = Surv(s.inv, delta) ~ mixture(age) + mixture(trt),
                              hazard = hazard.str,
                              hazardnodes = knot.seq[1:knot.k],
                              logscale = FALSE,
                              subject = 'ID',
                              data = YS,
                              ng = 3))

hare4 <- gridsearch(rep = 30, maxiter = 15, minit = hare1,
                    Jointlcmm(y ~ visit + age + sex + trt,
                              random =~ visit,
                              mixture =~ age + visit,
                              survival = Surv(s.inv, delta) ~ mixture(age) + mixture(trt),
                              hazard = hazard.str,
                              hazardnodes = knot.seq[1:knot.k],
                              logscale = FALSE,
                              subject = 'ID',
                              data = YS,
                              ng = 4))

#################
#Comparing model
#fit statistics
#################
summarytable(weibull1,hare1,weibull2,hare2,weibull3,hare3,weibull4,hare4)

###################
#Regression tables
###################
summary(weibull3)
summary(hare3)

#########
#N = 800
#########

X1B <- design_matrix(2704,400)
X2B <- design_matrix(2704,300)
X3B <- design_matrix(2704,100)

beta1B <- c(1.5, 0.5, 0.25)
beta2B <- c(0.5, -0.5, -0.25)
beta3B <- c(3.0, 1.0, 0.125)


gamma1B <- c(1.125, 0.25)
gamma2B <- c(-0.05, 0.55)
gamma3B <- c(0.35, -0.35)


#Scale
lambda1B <- 1
lambda2B <- 0.9
lambda3B <- 1.1

#Shape
psi1B <- 0.5
psi2B <- 1.5
psi3B <- 3


Y1B <- design_all(X1B, beta = beta1B, lambda = lambda1B, psi = psi1B, gamma = gamma1B, IDoffset = 0)
Y2B <- design_all(X2B, beta = beta2B, lambda = lambda2B, psi = psi2B, gamma = gamma2B, IDoffset = 400)
Y3B <- design_all(X3B, beta = beta3B, lambda = lambda3B, psi = psi3B, gamma = gamma3B, IDoffset = 700)

YSB <- data.frame(rbind(Y1B,Y2B,Y3B))

cov.list <- c("sex","age","trt")

hare.modB <- hare(data = YSB$s.inv, delta = YSB$delta, cov = YSB[,cov.list])
hare.modB$knots

knot.kB <- hare.modB$knots[1,1]
hazard.kB <- knot.kB + 2
hazard.strB <- paste0(hazard.kB,"-manual-splines")
knot.seqB <- c(hare.modB$knots[1,-1])

weibull1B <- Jointlcmm(y ~ visit + age + sex + trt,
                       random =~ visit,
                       survival = Surv(s.inv, delta) ~ age + sex + trt,
                       hazard = "Weibull",
                       logscale = FALSE,
                       subject = 'ID',
                       data = YSB,
                       ng = 1)

hare1B <- Jointlcmm(y ~ visit + age + sex + trt,
                    random =~ visit,
                    survival = Surv(s.inv, delta) ~ age + trt,
                    hazard = hazard.strB,
                    hazardnodes = knot.seqB[1:knot.kB],
                    logscale = FALSE,
                    subject = 'ID',
                    data = YSB,
                    ng = 1)

weibull2B <- Jointlcmm(y ~ visit + age + sex + trt,
                       random =~ visit,
                       mixture =~ age + visit,
                       survival = Surv(s.inv, delta) ~ mixture(age) + mixture(sex) + mixture(trt),
                       hazard = "Weibull",
                       logscale = FALSE,
                       subject = 'ID',
                       data = YSB,
                       B = weibull1B,
                       ng = 2)

weibull3B <- Jointlcmm(y ~ visit + age + sex + trt,
                       random =~ visit,
                       mixture =~ age + visit,
                       survival = Surv(s.inv, delta) ~ mixture(age) + mixture(sex) + mixture(trt),
                       hazard = "Weibull",
                       logscale = FALSE,
                       subject = 'ID',
                       data = YSB,
                       B = weibull1B,
                       ng = 3)

weibull4B <- Jointlcmm(y ~ visit + age + sex + trt,
                       random =~ visit,
                       mixture =~ age + visit,
                       survival = Surv(s.inv, delta) ~ mixture(age) + mixture(sex) + mixture(trt),
                       hazard = "Weibull",
                       logscale = FALSE,
                       subject = 'ID',
                       data = YSB,
                       B = weibull1B,
                       ng = 4)

hare2B <- gridsearch(rep = 30, maxiter = 15, minit = hare1B,
                     Jointlcmm(y ~ visit + age + sex + trt,
                               random =~ visit,
                               mixture =~ age + visit,
                               survival = Surv(s.inv, delta) ~ mixture(age) + mixture(trt),
                               hazard = hazard.strB,
                               hazardnodes = knot.seqB[1:knot.kB],
                               logscale = FALSE,
                               subject = 'ID',
                               data = YSB,
                               ng = 2))

hare3B <- gridsearch(rep = 50, maxiter = 20, minit = hare1B,
                     Jointlcmm(y ~ visit + age + sex + trt,
                               random =~ visit,
                               mixture =~ age + visit,
                               survival = Surv(s.inv, delta) ~ mixture(age) + mixture(trt),
                               hazard = hazard.strB,
                               hazardnodes = knot.seqB[1:knot.kB],
                               logscale = FALSE,
                               subject = 'ID',
                               data = YSB,
                               ng = 3))

hare4B <- gridsearch(rep = 30, maxiter = 15, minit = hare1B,
                     Jointlcmm(y ~ visit + age + sex + trt,
                               random =~ visit,
                               mixture =~ age + visit,
                               survival = Surv(s.inv, delta) ~ mixture(age) + mixture(trt),
                               hazard = hazard.strB,
                               hazardnodes = knot.seqB[1:knot.kB],
                               logscale = FALSE,
                               subject = 'ID',
                               data = YSB,
                               ng = 4))


summarytable(weibull1B,hare1B,weibull2B,hare2B,weibull3B,hare3B,weibull4B,hare4B)

summary(weibull3B)
summary(hare3B)

#########
#N = 1200
#########

X1C <- design_matrix(2702,600)
X2C <- design_matrix(2702,450)
X3C <- design_matrix(2702,150)

beta1C <- c(1.5, 0.5, 0.25)
beta2C <- c(0.5, -0.5, -0.25)
beta3C <- c(3.0, 1.0, 0.125)


gamma1C <- c(1.125, 0.25)
gamma2C <- c(-0.05, 0.55)
gamma3C <- c(0.35, -0.35)


#Scale
lambda1C <- 1
lambda2C <- 0.9
lambda3C <- 1.1

#Shape
psi1C <- 0.5
psi2C <- 1.5
psi3C <- 3


Y1C <- design_all(X1C, beta = beta1C, lambda = lambda1C, psi = psi1C, gamma = gamma1C, IDoffset = 0)
Y2C <- design_all(X2C, beta = beta2C, lambda = lambda2C, psi = psi2C, gamma = gamma2C, IDoffset = 600)
Y3C <- design_all(X3C, beta = beta3C, lambda = lambda3C, psi = psi3C, gamma = gamma3C, IDoffset = 1050)

YSC <- data.frame(rbind(Y1C,Y2C,Y3C))

cov.list <- c("sex","age","trt")

hare.modC <- hare(data = YSC$s.inv, delta = YSC$delta, cov = YSC[,cov.list])
hare.modC$knots

knot.kC <- hare.modC$knots[1,1]
hazard.kC <- knot.kC + 2
hazard.strC <- paste0(hazard.kC,"-manual-splines")
knot.seqC <- c(hare.modC$knots[1,-1])

weibull1C <- Jointlcmm(y ~ visit + age + sex + trt,
                       random =~ visit,
                       survival = Surv(s.inv, delta) ~ age + sex + trt,
                       hazard = "Weibull",
                       logscale = FALSE,
                       subject = 'ID',
                       data = YSC,
                       ng = 1)

hare1C <- Jointlcmm(y ~ visit + age + sex + trt,
                    random =~ visit,
                    survival = Surv(s.inv, delta) ~ age + trt,
                    hazard = hazard.strC,
                    hazardnodes = knot.seqC[1:knot.kC],
                    logscale = FALSE,
                    subject = 'ID',
                    data = YSC,
                    ng = 1)

weibull2C <- Jointlcmm(y ~ visit + age + sex + trt,
                       random =~ visit,
                       mixture =~ age + visit,
                       survival = Surv(s.inv, delta) ~ mixture(age) + mixture(sex) + mixture(trt),
                       hazard = "Weibull",
                       logscale = FALSE,
                       subject = 'ID',
                       data = YSC,
                       B = weibull1C,
                       ng = 2)

weibull3C <- Jointlcmm(y ~ visit + age + sex + trt,
                       random =~ visit,
                       mixture =~ age + visit,
                       survival = Surv(s.inv, delta) ~ mixture(age) + mixture(sex) + mixture(trt),
                       hazard = "Weibull",
                       logscale = FALSE,
                       subject = 'ID',
                       data = YSC,
                       B = weibull1C,
                       ng = 3)

weibull4C <- Jointlcmm(y ~ visit + age + sex + trt,
                       random =~ visit,
                       mixture =~ age + visit,
                       survival = Surv(s.inv, delta) ~ mixture(age) + mixture(sex) + mixture(trt),
                       hazard = "Weibull",
                       logscale = FALSE,
                       subject = 'ID',
                       data = YSC,
                       B = weibull1C,
                       ng = 4)

hare2C <- gridsearch(rep = 30, maxiter = 15, minit = hare1C,
                     Jointlcmm(y ~ visit + age + sex + trt,
                               random =~ visit,
                               mixture =~ age + visit,
                               survival = Surv(s.inv, delta) ~ mixture(age) + mixture(trt),
                               hazard = hazard.strC,
                               hazardnodes = knot.seqC[1:knot.kC],
                               logscale = FALSE,
                               subject = 'ID',
                               data = YSC,
                               ng = 2))

hare3C <- gridsearch(rep = 30, maxiter = 15, minit = hare1C,
                     Jointlcmm(y ~ visit + age + sex + trt,
                               random =~ visit,
                               mixture =~ age + visit,
                               survival = Surv(s.inv, delta) ~ mixture(age) + mixture(trt),
                               hazard = hazard.strC,
                               hazardnodes = knot.seqC[1:knot.kC],
                               logscale = FALSE,
                               subject = 'ID',
                               data = YSC,
                               ng = 3))

hare4C <- gridsearch(rep = 30, maxiter = 15, minit = hare1C,
                     Jointlcmm(y ~ visit + age + sex + trt,
                               random =~ visit,
                               mixture =~ age + visit,
                               survival = Surv(s.inv, delta) ~ mixture(age) + mixture(trt),
                               hazard = hazard.strC,
                               hazardnodes = knot.seqC[1:knot.kC],
                               logscale = FALSE,
                               subject = 'ID',
                               data = YSC,
                               ng = 4))


summarytable(weibull1C,hare1C,weibull2C,hare2C,weibull3C,hare3C,weibull4C,hare4C)

summary(weibull3C)
summary(hare3C)

######################
#Plots for manuscript
######################

## Figure 2

## By definition, three classes were produced as three separate datasets
## whose identification can be determined by ID

## Using ID to produce latent class groups 1 - 3
data.for.plotting <- YS %>% group_by(ID) %>%
 mutate(time = row_number(),
        real_lc = case_when(between(ID,1,200) ~ 1,
                            between(ID,201,350) ~ 2,
                            between(ID,351,400) ~ 3)) %>%
 ungroup()
 
lc_labels <- c(`1` = "True LC 1", `2` = "True LC 2", `3` = "True LC 3")
 
stimes <- data.for.plotting %>% ggplot(aes(x = s.inv)) +
 geom_histogram(color = 'black') + 
 theme_bw() +
 theme(panel.grid.minor = element_blank(),
       panel.grid.major = element_blank(),
       text = element_text(size = 20),
       legend.position = NULL) +
 xlab("Survival Times") +
 ylab("Count") +
 facet_wrap(. ~ real_lc,
            labeller = as_labeller(lc_labels))
 
longitudinal <- data.for.plotting %>% ggplot(aes(x = time, y = y, group = ID)) +
 geom_line() +
 theme_bw() +
 theme(panel.grid.minor = element_blank(),
       panel.grid.major = element_blank(),
       text = element_text(size = 20),
       legend.position = NULL) +
 xlab("Visit") +
 scale_y_continuous(limits = c(-2,6),
                    breaks = seq(-2, 6,0.5),
                    name = "Longitudinal Outcome") +
 facet_wrap(. ~ real_lc,
            labeller = as_labeller(lc_labels))

#Using patchwork to place plots into an array for publication.  
#If you want to inspect the separate plots, simply refer to `stimes` or `longitudinal` for the survival times, longitudinal trajectories, respectively

library(patchwork)
sd.histlines <- stimes / longitudinal
 

## Figure 3

#Extracting the hazard estimates from each three-class JLCMM
weib.surv <- weibull3$predSurv
hare.surv <- hare3$predSurv


#Calculating hazard function and survival probabilities (determined from cumulative hazard) for Weibull model
weib.surv.haz <- weib.surv %>% as.data.frame() %>% select(time,event1.RiskFct1,event1.RiskFct2,event1.RiskFct3) %>%
  pivot_longer(cols=c(event1.RiskFct1,event1.RiskFct2,event1.RiskFct3),names_to = "class",
               values_to = "hazard", names_pattern = "event1.RiskFct(.*)") %>% as.data.frame()

weib.surv.cumhaz <- weib.surv %>% as.data.frame() %>% select(time,event1.CumRiskFct1,event1.CumRiskFct2,event1.CumRiskFct3) %>%
  pivot_longer(cols=c(event1.CumRiskFct1,event1.CumRiskFct2,event1.CumRiskFct3),names_to = "class",
               values_to = "cumhazard", names_pattern = "event1.CumRiskFct(.*)") %>%
  mutate(surv = exp(-cumhazard)) %>% as.data.frame()

#Joining into one dataframe
weib.survival.pred.data <- left_join(weib.surv.haz,weib.surv.cumhaz,by=c("time","class"))

#Now calculating hazard function and survival probabilities for HARE model
hare.surv.haz <- hare.surv %>% as.data.frame() %>% select(time,event1.RiskFct1,event1.RiskFct2,event1.RiskFct3) %>%
  pivot_longer(cols=c(event1.RiskFct1,event1.RiskFct2,event1.RiskFct3),names_to = "class",
               values_to = "hazard", names_pattern = "event1.RiskFct(.*)") %>% as.data.frame()

hare.surv.cumhaz <- hare.surv %>% as.data.frame() %>% select(time,event1.CumRiskFct1,event1.CumRiskFct2,event1.CumRiskFct3) %>%
  pivot_longer(cols=c(event1.CumRiskFct1,event1.CumRiskFct2,event1.CumRiskFct3),names_to = "class",
               values_to = "cumhazard", names_pattern = "event1.CumRiskFct(.*)") %>%
  mutate(surv = exp(-cumhazard)) %>% as.data.frame()

#Joining into one dataframe
hare.survival.pred.data <- left_join(hare.surv.haz,hare.surv.cumhaz,by=c("time","class"))


#Predicting the longitudinal outcomes over each visit for both three-class JLCMMs
datnew <- data.frame(visit = seq(0,2,length = 100))
datnew$age <- 0
datnew$sex <- 0
datnew$trt <- 0
weib.predict <- predictY(weibull3, newdata = datnew, var.time = 'visit')
hare.predict <- predictY(hare3, newdata = datnew, var.time = 'visit')

weib.pred.data <- data.frame(weib.predict$times,weib.predict$pred) %>%
  pivot_longer(cols = c(Ypred_class1,Ypred_class2,Ypred_class3), names_to = "class",
               values_to = "predictY", names_pattern = "Ypred_class(.*)") %>% as.data.frame()

hare.pred.data <- data.frame(hare.predict$times,hare.predict$pred) %>%
  pivot_longer(cols = c(Ypred_class1,Ypred_class2,Ypred_class3), names_to = "class",
               values_to = "predictY", names_pattern = "Ypred_class(.*)") %>% as.data.frame()

## JLCMMs and other latent class models have label switching issues due to the labels being arbitrarily determined.
## Need to verify which assigned class refers to each latent class by using postprob(weibull3) and postprob(hare3).
## Can change assign_class variable as needed to refer to true latent class by determining whether postprob() classification overlaps with true classification at least 85%.
## In simulations, these overlaps were typically 95% or better.

## While machine for these simulations have consistently given the classifications in the code below, I have found on occassion that changing 
## computers has caused different label assignments (but not different proportions, posterior probabilities, or regression estimates)


#Longitudinal trajectories for each class
wb.prd <- weib.pred.data %>% mutate(assign_class = case_when(class == 1 ~ 2,
                                                             class == 2 ~ 1,
                                                             class == 3 ~ 3)) %>%
  ggplot(aes(x = visit, y = predictY)) + geom_line(aes(linetype = as.factor(assign_class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  scale_x_continuous(labels = c(1,2,3,4,5),
                     breaks = c(0,0.5,1,1.5,2),
                     name = "Visit") +
  scale_y_continuous(limits = c(-2,6),
                     breaks = seq(-2,6,0.5),
                     name = "Predicted Outcome") +
  labs(linetype = "Assigned Class")

hr.prd <- hare.pred.data %>% mutate(assign_class = case_when(class == 1 ~ 2,
                                                             class == 2 ~ 1,
                                                             class == 3 ~ 3)) %>%
  ggplot(aes(x = visit, y = predictY)) + geom_line(aes(linetype = as.factor(assign_class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  scale_x_continuous(labels = c(1,2,3,4,5),
                     breaks = c(0,0.5,1,1.5,2),
                     name = "Visit") +
  scale_y_continuous(limits = c(-2,6),
                     breaks = seq(-2,6,0.5),
                     name = "Predicted Outcome") +
  labs(linetype = "Assigned Class")

#Weibull hazard functions for each class
wb.haz <- weib.survival.pred.data %>% mutate(assign_class = case_when(class == 1 ~ 2,
                                                                      class == 2 ~ 1,
                                                                      class == 3 ~ 3)) %>%
  ggplot(aes(x = time, y = hazard)) + geom_path(aes(linetype = as.factor(assign_class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Time") +
  scale_y_continuous(limits = c(0,20),
                     breaks = seq(0,20,2.0),
                     expand = c(0,0),
                     name = "Predicted Hazard") +
  labs(linetype = "Assigned Class")

#Weibull survival probabilities for each class
wb.surv <- weib.survival.pred.data %>% mutate(assign_class = case_when(class == 1 ~ 2,
                                                                       class == 2 ~ 1,
                                                                       class == 3 ~ 3)) %>%
  ggplot(aes(x = time, y = surv)) + geom_path(aes(linetype = as.factor(assign_class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Time") +
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0),
                     name = "Predicted Survival Probability") +
  labs(linetype = "Assigned Class")

#HARE hazard functions for each class
hr.haz <- hare.survival.pred.data %>% mutate(assign_class = case_when(class == 1 ~ 2,
                                                                      class == 2 ~ 1,
                                                                      class == 3 ~ 3)) %>%
  ggplot(aes(x = time, y = hazard)) + geom_path(aes(linetype = as.factor(assign_class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Time") +
  scale_y_continuous(limits = c(0,140),
                     breaks = seq(0,140,20.0),
                     expand = c(0,0),
                     name = "Predicted Hazard") +
  labs(linetype = "Assigned Class")

#HARE survival probabilities for each class
hr.surv <- hare.survival.pred.data %>% mutate(assign_class = case_when(class == 1 ~ 2,
                                                                       class == 2 ~ 1,
                                                                       class == 3 ~ 3)) %>%
  ggplot(aes(x = time, y = surv)) + geom_path(aes(linetype = as.factor(assign_class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Time") +
  scale_y_continuous(limits = c(0,1,0.1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0),
                     name = "Predicted Survival Probability") +
  labs(linetype = "Assigned Class")

#Patchwork again
wb_hr_preds <- (wb.prd + wb.haz + wb.surv) / (hr.prd + hr.haz + hr.surv) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
