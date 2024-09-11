#Note: Different Versions of 'lcmm' and R may affect the results.  Refer to the manuscript for appropriate versions to replicate analyses done here.
#Note: Code here developed for original submission of paper.  The code itself has not changed between versions, but some tables are now removed
## from the manuscript and left in the online companion (09-01-2024).
library(MASS)
library(tidyverse)
library(ggplot2)
library(knitr)
library(kableExtra)
library(dplyr)
library(polspline)
library(survival)
library(magrittr)
library(combinat)
library(lcmm)
library(NormPsy)

set.seed(2702)

paquid <- lcmm::paquid #Longitudinal data on cognitive and physical againg in the elderly


paquid$normMMSE <- NormPsy::normMMSE(paquid$MMSE)
paquid$age65 <- (paquid$age-65)/10
head(paquid) #Check.

paquidS <- paquid[paquid$agedem > paquid$age_init, ]
paquidS <- na.omit(paquidS) #Subset nonmissing data for hare() to run


paquidS.2 <- paquidS %>% mutate(stime = agedem - age_init) %>% select(ID,stime,dem,CEP,male) %>%
  group_by(ID) %>% unique() %>% ungroup() #Using right-censored time-to-event dataframe for hare() to run

cov.list <- c("CEP","male")

data.for.hare <- paquidS.2 %>% as.data.frame() #In data.frame form to be sed to estimate knots. Jointlcmm() use paquidS data.

hare.paquid <- hare(data=data.for.hare$stime, delta = data.for.hare$dem, cov = data.for.hare[,cov.list])

hare.paquid$knots #18.82132, 19.15414, 19.95; add 65 to get hazard nodes in paquidS

mj1 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 random =~ poly(age65, degree = 2, raw = TRUE),
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "Weibull", subject = 'ID', data = paquidS, ng = 1)
mj2 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 mixture =~ poly(age65, degree = 2, raw = TRUE),
                 random =~ poly(age65, degree = 2, raw = TRUE),
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "Weibull", subject = 'ID', data = paquidS, ng = 2, B = mj1)
mj3 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 mixture =~ poly(age65, degree = 2, raw = TRUE),
                 random =~ poly(age65, degree = 2, raw = TRUE),
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "Weibull", subject = 'ID', data = paquidS, ng = 3, B = mj1)
mj4 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 mixture =~ poly(age65, degree = 2, raw = TRUE),
                 random =~ poly(age65, degree = 2, raw = TRUE),
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "Weibull", subject = 'ID', data = paquidS, ng = 4, B = mj1)

hj1 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 random =~ poly(age65, degree = 2, raw = TRUE),
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "5-manual-splines", hazardnodes = c(83.82132, 84.15414, 84.95), logscale = TRUE,
                 subject = 'ID', data = paquidS, ng = 1) 

hj2 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 mixture =~ poly(age65, degree = 2, raw = TRUE),
                 random =~ poly(age65, degree = 2, raw = TRUE),
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "5-manual-splines", hazardnodes = c(83.82132, 84.15414, 84.95), hazardtype = "PH", logscale = TRUE,
                 subject = 'ID', data = paquidS, ng = 2, B = hj1)
hj3 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 mixture =~ poly(age65, degree = 2, raw = TRUE),
                 random =~ poly(age65, degree = 2, raw = TRUE),
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "5-manual-splines", hazardnodes = c(83.82132, 84.15414, 84.95), hazardtype = "PH", logscale = TRUE,
                 subject = 'ID', data = paquidS, ng = 3, B = hj1)
hj4 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 mixture =~ poly(age65, degree = 2, raw = TRUE),
                 random =~ poly(age65, degree = 2, raw = TRUE),
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "5-manual-splines", hazardnodes = c(83.82132, 84.15414, 84.95), hazardtype = "PH", logscale = TRUE,
                 subject = 'ID', data = paquidS, ng = 4, B = hj1)


## Among these models, three-class HARE solution has lowest BIC (15536.37). 
#
## Further, Weibull-determined solutions have some difficulty finding
##  solutions to the three- and four-class solutions.  HARE also has
##  some difficulty converging to four-class solution.
summarytable(mj1, hj1, mj2, hj2, mj3, hj3, mj4, hj4)



##########################################
#Using manual starting values
# and gridsearch() to help
# find global maxima.  Weibull
# solutions used starting values
# based on Proust-Lima et al. (2017).
#
#While some manual starting values
# were used for HARE models, there
# was some difficulty in convergence.
# Therefore, gridsearch() was used
# instead.
##########################################

Binit <- rep(0, length(mj2$best) + 6)
Binit[c(2, 5:10, 12, 13, 15, 16, 18, 19:(length(Binit)))] <- mj2$best #Using best values from two-class solution as starting values

Binit[c(1, 3, 4, 11, 14, 17)] <- c(0, 0.11, 4, 70, 0, 0) #Changing a few of these starting values to new ones.  Also need to accommodate additional parameters for three classes

mj3b <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                  mixture =~ poly(age65, degree = 2, raw = TRUE),
                  random =~ poly(age65, degree = 2, raw = TRUE),
                  survival = Surv(age_init, agedem, dem) ~ CEP + male,
                  hazard = "Weibull", subject = 'ID', data = paquidS, ng = 3, B = Binit)


Binit.2 <- rep(0, length(mj3b$best)+6)
Binit.2[c(1, 2, 4:7, 10:15, 17:19, 21:23, 25:length(Binit.2))] <- mj3b$best
Binit.2[c(3, 8, 9, 16, 20, 24)] <- c(0, 0.1, 10, 60, 5, -10)

mj4b <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                  mixture =~ poly(age65, degree = 2, raw = TRUE),
                  random =~ poly(age65, degree = 2, raw = TRUE),
                  survival = Surv(age_init, agedem, dem) ~ CEP + male,
                  hazard = "Weibull", subject = 'ID', data = paquidS, ng = 4, B = Binit.2)

#Grid search for HARE models
hj3c <- gridsearch(rep = 30, maxiter = 15, minit = hj1,
                   Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                             mixture =~ poly(age65, degree = 2, raw = TRUE),
                             random =~ poly(age65, degree = 2, raw = TRUE),
                             survival = Surv(age_init, agedem, dem) ~ CEP + male,
                             hazard = "5-manual-splines", hazardnodes = c(83.82132, 84.15414, 84.95), hazardtype = "PH", logscale = TRUE,
                             subject = 'ID', data = paquidS, ng = 3))

hj4c <- gridsearch(rep = 30, maxiter = 15, minit = hj1,
                   Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                             mixture =~ poly(age65, degree = 2, raw = TRUE),
                             random =~ poly(age65, degree = 2, raw = TRUE),
                             survival = Surv(age_init, agedem, dem) ~ CEP + male,
                             hazard = "5-manual-splines", hazardnodes = c(83.82132, 84.15414, 84.95), hazardtype = "PH", logscale = TRUE,
                             subject = 'ID', data = paquidS, ng = 4))

#Grid search with class-specific hazard types for HARE
hj2sp <- gridsearch(rep = 30, maxiter = 15, minit = hj1,
                    Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                              mixture =~ poly(age65, degree = 2, raw = TRUE),
                              random =~ poly(age65, degree = 2, raw = TRUE),
                              survival = Surv(age_init, agedem, dem) ~ CEP + male,
                              hazard = "5-manual-splines", hazardnodes = c(83.82132, 84.15414, 84.95), logscale = TRUE,
                              subject = 'ID', data = paquidS, ng = 2)) 

hj3sp <- gridsearch(rep = 30, maxiter = 15, minit = hj1,
                    Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                              mixture =~ poly(age65, degree = 2, raw = TRUE),
                              random =~ poly(age65, degree = 2, raw = TRUE),
                              survival = Surv(age_init, agedem, dem) ~ CEP + male,
                              hazard = "5-manual-splines", hazardnodes = c(83.82132, 84.15414, 84.95), logscale = TRUE,
                              subject = 'ID', data = paquidS, ng = 3))

hj4sp <- gridsearch(rep = 30, maxiter = 15, minit = hj1,
                    Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                              mixture =~ poly(age65, degree = 2, raw = TRUE),
                              random =~ poly(age65, degree = 2, raw = TRUE),
                              survival = Surv(age_init, agedem, dem) ~ CEP + male,
                              hazard = "5-manual-splines", hazardnodes = c(83.82132, 84.15414, 84.95), logscale = TRUE,
                              subject = 'ID', data = paquidS, ng = 4))


summarytable(mj1, hj1, mj2, hj2, hj2sp, mj3b, hj3, hj3sp, mj4b, hj4, hj4sp) #Table 8

postprob(mj4b) #Table 9
summary(mj4b) #Table 13

postprob(hj3) #Table 10
summary(hj3) #Table 14



####################
#Empircal Functions
# (Figure 4)
####################

### Setting up data for empirical hazard and empirical survival curves
for.emp.ests <- paquidS %>% mutate(stime = agedem - age_init) %>% select(ID,age_init,stime,dem,CEP,male) %>%
  group_by(ID) %>% unique() %>% ungroup() %>% as.data.frame()

#Values below recreate Table 7

N <- numeric()
c <- numeric()
h <- numeric()
S <- numeric()
for (t in seq(0,max(for.emp.ests$stime),2)){
  nn <- for.emp.ests %>% group_by(stime >= t) %>% summarise(freq = n()) %>% filter(.[,1] == TRUE) %>% pull()
  cc <- for.emp.ests %>% group_by(stime >= t & stime < t + 2 & dem == 1) %>% summarise(freq = n()) %>% filter(.[,1] == TRUE) %>% pull() #Used two-year interval since sufficiently small
  hh <- cc/nn
  ss <- nn/490
  N <- c(N,nn)
  c <- c(c,cc)
  h <- c(h,hh) #empirical hazard
  S <- c(S,ss) #empirical survival function
}
c <- c(c,0) #padding 0 for last interval
h <- c(h,0) #padding 0 for last interval
N #People not censored or experienced event by end of each interval
c #Number of deaths in interval
h #Empirical hazard
S #Empirical Survival
t <- seq(0,max(for.emp.ests$stime),2)

data.for.emp <- data.frame(t,h,S)

#Empirical hazard plot
aq <- data.for.emp %>% ggplot(aes(x = t, y = h)) + geom_line() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Age in Years Since 65") +
  scale_y_continuous(name = "Empirical Hazard")

#Empirical Survival Probability plot
bq <- data.for.emp %>% ggplot(aes(x = t, y = S)) + geom_line() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Age in Years Since 65") +
  scale_y_continuous(name = "Empirical Survival Probability")

library(patchwork) #Using patchwork here to arrange the plots
qq <- aq / bq

####################################
#Comparing latent class assignment
# (Table 11)
####################################

hj3c.pprobs <- hj3$pprob

mj4b.pprobs <- mj4b$pprob

names(hj3c.pprobs) <- c("ID","class_hare","pprob1_hare","pprob2_hare","pprob3_hare")
names(mj4b.pprobs) <- c("ID","class_weib","pprob1_weib","pprob2_weib","pprob3_weib","pprob4_weib")

all_classes <- left_join(hj3c.pprobs,mj4b.pprobs,by="ID")

all_classes %$% table(class_hare,class_weib)


########################
#Creating Figure 5
########################
mj4b.surv <- mj4b$predSurv
hj3c.surv <- hj3$predSurv

mj4b.surv.haz <- mj4b.surv %>% as.data.frame() %>% select(time,event1.RiskFct1,event1.RiskFct2,event1.RiskFct3,event1.RiskFct4) %>%
  pivot_longer(cols=c(event1.RiskFct1,event1.RiskFct2,event1.RiskFct3,event1.RiskFct4), names_to = "class",
               values_to = "hazard", names_pattern = "event1.RiskFct(.*)") %>% as.data.frame()

mj4b.surv.cumhaz <- mj4b.surv %>% as.data.frame() %>% select(time,event1.CumRiskFct1,event1.CumRiskFct2,event1.CumRiskFct3,event1.CumRiskFct4) %>%
  pivot_longer(cols=c(event1.CumRiskFct1,event1.CumRiskFct2,event1.CumRiskFct3,event1.CumRiskFct4), names_to = "class",
               values_to = "cumhazard", names_pattern = "event1.CumRiskFct(.*)") %>% 
  mutate(surv = exp(-cumhazard)) %>% as.data.frame()

hj3c.surv.haz <- hj3c.surv %>% as.data.frame() %>% select(time,event1.RiskFct1,event1.RiskFct2,event1.RiskFct3) %>%
  pivot_longer(cols=c(event1.RiskFct1,event1.RiskFct2,event1.RiskFct3), names_to = "class",
               values_to = "hazard", names_pattern = "event1.RiskFct(.*)") %>% as.data.frame()

hj3c.surv.cumhaz <- hj3c.surv %>% as.data.frame() %>% select(time,event1.CumRiskFct1,event1.CumRiskFct2,event1.CumRiskFct3) %>%
  pivot_longer(cols=c(event1.CumRiskFct1,event1.CumRiskFct2,event1.CumRiskFct3), names_to = "class",
               values_to = "cumhazard", names_pattern = "event1.CumRiskFct(.*)") %>% 
  mutate(surv = exp(-cumhazard)) %>% as.data.frame()

mj4b.survival.pred.data <- left_join(mj4b.surv.haz,mj4b.surv.cumhaz,by=c("time","class"))
hj3c.survival.pred.data <- left_join(hj3c.surv.haz,hj3c.surv.cumhaz,by=c("time","class"))

datnew <- data.frame(age65 = seq(0, 3, length=100))
datnew$male <- 0
datnew$CEP <- 0
par(mfrow = c(1, 2))
mj4b.pred <- predictY(mj4b, newdata = datnew, var.time = "age65")
hj3c.pred <- predictY(hj3, newdata = datnew, var.time = "age65") 

mj4b.pred.data <- data.frame(mj4b.pred$times, mj4b.pred$pred) %>%
  pivot_longer(cols = c(Ypred_class1, Ypred_class2, Ypred_class3, Ypred_class4), names_to = "class",
               values_to = "predictY", names_pattern = "Ypred_class(.*)") %>% as.data.frame()

hj3c.pred.data <- data.frame(hj3c.pred$times, hj3c.pred$pred) %>%
  pivot_longer(cols = c(Ypred_class1, Ypred_class2, Ypred_class3), names_to = "class",
               values_to = "predictY", names_pattern = "Ypred_class(.*)") %>% as.data.frame()

mj4b.pred.plot <- mj4b.pred.data %>% ggplot(aes(x = age65, y = predictY)) + 
  geom_line(aes(linetype = as.factor(class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Decades from Age of 65") +
  ylab("Normed MMSE Score") +
  # scale_color_manual(values = c("black", "red4", "yellowgreen", "steelblue4")) +
  scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotted")) +
  labs(color = "Latent Class", linetype = "Latent Class") +
  annotate("text",x = 0, y = 83, label = "(A)", size = 6)

mj4b.surv.plot <- mj4b.survival.pred.data %>% ggplot(aes(x = time, y = surv)) + 
  geom_line(aes(linetype = as.factor(class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Age (years)") +
  ylab("Dementia-Free Probability") +
  # scale_color_manual(values = c("black", "red4", "yellowgreen", "steelblue4")) +
  scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotted")) +
  labs(color = "Latent Class", linetype = "Latent Class") +
  annotate("text",x = 65, y = 1.1, label = "(B)", size = 6)

hj3c.pred.plot <- hj3c.pred.data %>% ggplot(aes(x = age65, y = predictY)) + 
  geom_line(aes(linetype = as.factor(class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Decades from Age of 65") +
  ylab("Normed MMSE Score") +
  # scale_color_manual(values = c("black", "red4", "yellowgreen")) +
  scale_linetype_manual(values = c("solid", "dashed", "twodash")) +
  labs(color = "Latent Class", linetype = "Latent Class") +
  annotate("text",x = 0, y = 83, label = "(C)", size = 6)

hj3c.surv.plot <- hj3c.survival.pred.data %>% ggplot(aes(x = time, y = surv)) + 
  geom_line(aes(linetype = as.factor(class)), size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        legend.position = 'bottom') +
  xlab("Age (years)") +
  ylab("Dementia-Free Probability") +
  # scale_color_manual(values = c("black", "red4", "yellowgreen")) +
  scale_linetype_manual(values = c("solid", "dashed", "twodash")) +
  labs(color = "Latent Class", linetype = "Latent Class") +
  annotate("text",x = 65, y = 1.1, label = "(D)", size = 6)

gg <- (mj4b.pred.plot + mj4b.surv.plot) / (hj3c.pred.plot + hj3c.surv.plot)



