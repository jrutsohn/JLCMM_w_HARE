library(dplyr)
library(ggplot2)
library(magrittr)
library(lcmm)
library(polspline)
library(survival)
library(mvtnorm)

`%ni%` <- Negate(`%in%`)

set.seed(2702)

your_directory <- setwd('') #SET THE DIRECTORY WHERE YOU WANT TO STORE THINGS HERE.  MAY NEED TO CREATE SUBFOLDERS.


#Assuming N = 400 for sample.  As seen in supplemental tables, increasing to N = 800 or N = 1200 while fixing latent class proportions 
## did not change parameter estimation, and standard error decreased.  To prevent disproportionately large or small LCs, the first LC
## proportion is sampled between 0.1 & 0.9.  The 1000 replicates of probabilities are generated first for inspection and later 
## evaluation.  This method randomly samples two two-digit numbers between 0.1 and 0.9, sorts them, and then calculates differences
## between the two numbers as well as the second number from 1.  This produces three probabilities for the multinomial sampling
## of the three latent classes.  The code below is adapted from a clever solution by Allan Cameron at stackoverflow here:
## https://stackoverflow.com/questions/72657010/generating-randomly-the-probabilities-that-sum-to-1 (Retrieved 07-07-24).



lc.probabilities <- t(replicate(1000, diff(c(0, sort(sample(seq(0.1, 0.9, 0.025), 2)), 1))))
lc.classes <- matrix(nrow = dim(lc.probabilities)[1], ncol = dim(lc.probabilities)[2], byrow = TRUE)
for(i in 1:dim(lc.probabilities)[1]){
  lc.classes[i,] <- rmultinom(n = 1, size = 400, prob = lc.probabilities[i,])
}


#As with the probabilities, the parameters for each model are randomly generated first for inspection and evaluation.  This should
## also decrease computation time as it minimizes for-loop processes.  Based on comments from reviewers, the simulated model has been
## simplified.  Various qualities of HARE to handle variables within the JLCMM can be seen viewing simulated results for the N = 400, 800, 1200 cases.
## In these simulations, we have three variables: visit (with random intercept), sex (with class-specific effects), and age (with class-specific effects
## in survival model and sample-wide effects in longitudinal model).  Class-specific intercepts (with person-level random intercept) are also included in the model.
## No variables have '0' effect in these simulations.  Weibull parameters are also randomly generated for each LC.  Beta are the mixed-model coefficients,
## theta are the survival coefficients, omega is Weibull shape, and phi is Weibull scale.

#MM Intercept
beta_0 <- t(replicate(1000,round(runif(3,-3,3),digits = 4)))
#MM Visit
beta_1 <- t(replicate(1000,round(runif(3,-3,3),digits = 4)))
#MM Sex
beta_2 <- t(replicate(1000,round(runif(3,-3,3),digits = 4)))
#MM Age
beta_3 <- replicate(1000,round(runif(1,-3,3),digits=4)) #Sampling one per replicate since the effect is sample-wide rather than class-specific

#Survival Sex
theta_1 <- t(replicate(1000,round(runif(3,-3,3),digits = 4)))
#Survival Age
theta_2 <- t(replicate(1000,round(runif(3,-3,3),digits = 4)))

#Weibull Shape
omega <- t(replicate(1000,round(runif(3,0,3),digits=4))) #Strictly positive
#Weibull Scale
phi <- t(replicate(1000,round(runif(3,0.5,1.5),digits=4))) #Strictly positive; placed in reasonable scale interval--too easy to distinguish LCs if scale is wildly different here.  This is true because the scale is a denominator value of the Weibull PDF.

#Need to create the design matrix for sex & age for one replicate to incorporate in for loop over all replicates
design_matrix <- function(n){
  sex <- ifelse(runif(n) < 0.5, 0, 1)
  #Allows for proportions of age distributions to be integers
  age_n1 <- round(n*0.2)
  age_n2 <- 2*age_n1
  age_n3 <- n - (age_n1+age_n2) #Remainder
  age_setup <- c(runif(age_n1, min = 30, max = 65), runif(age_n2, min = 65, max = 75), runif(age_n3, min = 75, max = 85))
  age <- (age_setup - 70)/10 #Centered & Scaled
  X <- cbind(sex,age)
  colnames(X) <- c("sex","age")
  return(X)
}

all_design_matrices <- list()
for(i in 1:1000){
  out <- c()
  vbeta_0 <- beta_0[i,]; vbeta_1 <- beta_1[i,]; vbeta_2 <- beta_2[i,]; vbeta_3 <- beta_3[i];
  vtheta_1 <- theta_1[i,]; vtheta_2 <- theta_2[i,]; 
  vomega <- omega[i,]; vphi <- phi[i,];
  for(j in 1:3){
    LC <- j
    beta0 <- vbeta_0[j]; beta1 <- vbeta_1[j]; beta2 <- vbeta_2[j]; beta3 <- vbeta_3;
    theta1 <- vtheta_1[j]; theta2 <- vtheta_2[j];
    omega_d <- vomega[j]; phi_d <- vphi[j];
    out[[j]] <- cbind(design_matrix(lc.classes[i,j]),LC,beta0,beta1,beta2,beta3,theta1,theta2,omega_d,phi_d)
  }
   all_design_matrices[[i]] <- data.frame(rbind(out[[1]],out[[2]],out[[3]])) %>% 
     mutate(ID = row_number()) %>% select(ID,everything())
}


#Creating both survival times and longitudinal times with this function.
## The idea is that the number of visits will be variable by individual (between 1 - 5),
## and this should reflect the survival times as well.  Survival times simulated with
## cumulative hazard inversion method.

#For the random effects in the mixed-effects model
sigma.2 <- rbind(c(0.04^2, 0.02*0.04^2), c(0.02*0.04^2, 0.02^2))

survival_times <- function(x){
  u <- runif(dim(x)[1],0,1)
  bi <- data.frame(rmvnorm(dim(x)[1], mean = rep(0, 2), sigma = sigma.2))
  names(bi) <- c("bi1","bi2")
  times <- sample(1:5,dim(x)[1],replace = TRUE)
  times.out <- data.frame(x,times,bi)
  total.sessions <- data.frame(cbind(sort(rep(seq(1,400,1),5)),rep(seq(1,5),400)))
  names(total.sessions) <- c("ID","session")
  times.out2 <- left_join(total.sessions,times.out,by="ID") %>% 
    group_by(ID) %>% 
    mutate(visit_flag = case_when(session <= times ~ 0,
                                  session > times ~ 1)) %>%
    filter(visit_flag == 0) %>%
    ungroup() %>%
    mutate(visit = 0.5*session - 0.5, #Putting on scale from 0 to 2.
           mu = (bi1 + beta0) +  (bi2 + beta1)*visit + beta2*sex + beta3*age,
           y = mu + rnorm(n(),0,0.1)) %>% #Mean + epsilon
    select(-times,-visit_flag)
  
  surv <- times.out2 %>% group_by(ID) %>% filter(row_number() == n()) %>% ungroup() %>%
    mutate(xtb = theta1*sex + theta2*age + 0.05*mu, #0.05 association coefficient
           sinv = (-log(u)/(phi_d*exp(xtb)))^(1/omega_d),
           delta.v = I(sinv < 2.1),
           s.inv = ifelse(delta.v,sinv,2.1),
           delta = as.numeric(delta.v)) %>%
    select(ID,s.inv,delta)
  
  jlcmm_data <-left_join(times.out2,surv,by="ID")
  
  
  return(jlcmm_data)
}


## Creating all the design matrices
jlcmm.data.list <- lapply(all_design_matrices,survival_times)

## Assigning a number to each replicate (this simplifies a lot of debugging)
for(i in 1:length(jlcmm.data.list)){
  jlcmm.data.list[[i]] %<>% mutate(simul_num = i)
}


#Covariate list for HARE
cov.list <- c("sex","age")

jlcmm_comparisons <- function(x){
  
  data.use <- data.frame(x)
  simul_num <- data.use %>% filter(row_number() == 1) %>% select(simul_num) %>% pull() #Simulation number.  Used for various row operations and output.
  
  #Running the HARE model to determine the hazard function and its knots
  hare.mod <- hare(data = data.use$s.inv, delta = data.use$delta, cov = data.use[,cov.list]) #s.inv is the time-to-event, delta is the event indicator (1 = experienced, 0 = censored)
  
  #Hazard regression knot extraction
  knot.k <- hare.mod$knots[1,1]
  hazard.k <- knot.k + 2
  hazard.str <- paste0(hazard.k,"-manual-splines") #For the JLCMM function
  knot.seq <- c(hare.mod$knots[1,-1])
  
  #One-class solutions for both Weibull and HARE models.  Need for comparisons with multiple-class solutions (and using parameter solutions from 1-class Weibull to initialize parameter estimates in multi-class solutions) 
  weibull1 <- Jointlcmm(y ~ visit + sex + age,
                        random =~ 1 + visit,
                        survival = Surv(s.inv, delta) ~ sex + age,
                        hazard = "Weibull",
                        logscale = FALSE,
                        subject = 'ID',
                        data = data.use,
                        ng = 1)
  
  hare1 <- Jointlcmm(y ~ visit + sex + age,
                     random =~ 1 + visit,
                     survival = Surv(s.inv, delta) ~ sex + age,
                     hazard = hazard.str,
                     hazardnodes = knot.seq[1:knot.k],
                     logscale = FALSE,
                     subject = 'ID',
                     data = data.use,
                     ng = 1)
  
  #Two-class solutions
  
  weibull2 <- Jointlcmm(y ~ visit + sex + age,
                        random =~ 1 + visit,
                        mixture =~ visit + sex,
                        survival = Surv(s.inv, delta) ~ mixture(sex) + mixture(age),
                        hazard = "Weibull",
                        logscale = FALSE,
                        subject = 'ID',
                        data = data.use,
                        B = weibull1,
                        ng = 2)
  
  hare2 <- gridsearch(rep = 30, maxiter = 15, minit = hare1,
                      Jointlcmm(y ~ visit + sex + age,
                                 random =~ 1 + visit,
                                 mixture =~ visit + sex,
                                 survival = Surv(s.inv, delta) ~ mixture(sex) + mixture(age),
                                 hazard = hazard.str,
                                 hazardnodes = knot.seq[1:knot.k],
                                 logscale = FALSE,
                                 subject = 'ID',
                                 data = data.use,
                                 ng = 2))
  
  #Three-class solutions
  
  weibull3 <- Jointlcmm(y ~ visit + sex + age,
                        random =~ 1 + visit,
                        mixture =~ visit + sex,
                        survival = Surv(s.inv, delta) ~ mixture(sex) + mixture(age),
                        hazard = "Weibull",
                        logscale = FALSE,
                        subject = 'ID',
                        data = data.use,
                        B = weibull1,
                        ng = 3)
  
  hare3 <- gridsearch(rep = 30, maxiter = 15, minit = hare1,
                      Jointlcmm(y ~ visit + sex + age,
                                random =~ 1 + visit,
                                mixture =~ visit + sex,
                                survival = Surv(s.inv, delta) ~ mixture(sex) + mixture(age),
                                hazard = hazard.str,
                                hazardnodes = knot.seq[1:knot.k],
                                logscale = FALSE,
                                subject = 'ID',
                                data = data.use,
                                ng = 3))
  
  #Four-class solutions: Can go higher if desired, but will abstain for the purposes of computation time
  
  weibull4 <- Jointlcmm(y ~ visit + sex + age,
                        random =~ 1 + visit,
                        mixture =~ visit + sex,
                        survival = Surv(s.inv, delta) ~ mixture(sex) + mixture(age),
                        hazard = "Weibull",
                        logscale = FALSE,
                        subject = 'ID',
                        data = data.use,
                        B = weibull1,
                        ng = 4)
  
  hare4 <- gridsearch(rep = 30, maxiter = 15, minit = hare1,
                      Jointlcmm(y ~ visit + sex + age,
                                random =~ 1 + visit,
                                mixture =~ visit + sex,
                                survival = Surv(s.inv, delta) ~ mixture(sex) + mixture(age),
                                hazard = hazard.str,
                                hazardnodes = knot.seq[1:knot.k],
                                logscale = FALSE,
                                subject = 'ID',
                                data = data.use,
                                ng = 4))
  
  #Sorting the class sizes to determine how to match classes.  Label switching inherently occurs since the labels are arbitrary, so
  # a way to match the assigned latent class with the true one was determined.  In the cases below, only the 3-class solution was
  # needed since true number of latent classes is alway 3.  If the best model is not the three class solution, then performance
  # measures are listed as missing anyway.
  
  sizes.sorted <- sort(lc.classes[simul_num,],index.return=TRUE,decreasing=TRUE)$x
  classes.sorted <- sort(lc.classes[simul_num,],index.return=TRUE,decreasing=TRUE)$ix
  
  weibull.sorted <- sort(table(weibull3$pprob$class),decreasing=TRUE)
  hare.sorted <- sort(table(hare3$pprob$class),decreasing=TRUE)
  
  weibull.truelc <- data.frame(classes.sorted,names(weibull.sorted))
  names(weibull.truelc) <- c("True_LC","Assigned_LC")
  hare.truelc <- data.frame(classes.sorted,names(hare.sorted))
  names(hare.truelc) <- c("True_LC","Assigned_LC")
  
  #Sort by Assigned_LC (ascending) then use True_LC for the replacement label.
  
  weibull.truelc %<>% arrange(Assigned_LC,ascending=TRUE)
  hare.truelc %<>% arrange(Assigned_LC,ascending=TRUE)
  
  #For label assignment
  weibull.pattern <- c(rep(weibull.truelc$True_LC,5),rep("",5))
  hare.pattern <- c(rep(hare.truelc$True_LC,5),rep("",5))
  
  #Creating the model selection summary table from which the best comparative models are decided (lowest BIC within Weibull/HARE models determines best fitting models)
  sum.table <- summarytable(weibull1,hare1,weibull2,hare2,weibull3,hare3,weibull4,hare4)
  best.weibull.model <- names(which.min(sum.table[c(grep("weibull",names(sum.table[,4]))),4]))
  best.hare.model <- names(which.min(sum.table[c(grep("hare",names(sum.table[,4]))),4]))
  
  #True parameters
  real.params <- matrix(c(beta_3[simul_num],t(beta_0[simul_num,]),t(beta_2[simul_num,]),t(beta_1[simul_num,]),t(theta_2[simul_num,]),t(theta_1[simul_num,])),ncol=1)

  
  #Total knots in HARE models are (hazard.k+2)*LC + LC-1  Therefore, need to extract best[-c(1:(hazard.k+2)*3 + 2)] from best HARE model (when LC = 3).
  ## For Weibull, will be best[-c(1:8)] no matter what (since all true number of LCs = 3).
  
  #Code below is a bit ad hoc, but the idea is that given non-convergence of the Hessian matrix (and nonestimability of SE), the models will be tested again
  ## after trying a broader grid search.  This can be done numerous times, but I have left it at two passes for both models just for simplicity and time.
  if(best.weibull.model == 'weibull3'){
    
    weibull.correct.model <- TRUE
    weibull.se <- matrix(sqrt(diag(VarCov(eval(parse(text=best.weibull.model)))))[-c(1:8)])
    if(!is.na(weibull.se[1])){
      weibull.best <- data.frame(names(eval(parse(text=best.weibull.model))$best[-c(1:8)]),eval(parse(text=best.weibull.model))$best[-c(1:8)],weibull.pattern,c(rep("surv",6),rep("long",10),rep("",4)),weibull.se)
      names(weibull.best) <- c("var","val","LC","model","se")
      weibull.best %<>% mutate(param = gsub("class\\d","",var),newvar = paste0(model," ",gsub("class\\d",paste0("class"),var),LC)) %>% arrange(model,param,LC,ascending = TRUE) %>% filter(var %ni% c('stderr','varcov 1', 'varcov 2', 'varcov 3'))
      
      weibull.comparison.data <- data.frame(weibull.best,real.params) %>% mutate(bias = val - real.params,
                                                                                 coverage = I(abs(val - real.params) < 2*se))
      
      weibull.se.estimated <- TRUE
    }
    if(is.na(weibull.se[1])){
      weibull3b <- gridsearch(rep = 50, maxiter = 25, minit = weibull1,
                          Jointlcmm(y ~ visit + sex + age,
                            random =~ 1 + visit,
                            mixture =~ visit + sex,
                            survival = Surv(s.inv, delta) ~ mixture(sex) + mixture(age),
                            hazard = "Weibull",
                            logscale = FALSE,
                            subject = 'ID',
                            data = data.use,
                            ng = 3))
      weibull2.se <- matrix(sqrt(diag(VarCov(weibull3b)))[-c(1:8)])
      if(!is.na(weibull2.se[1])){
        weibull.best <- data.frame(names(weibull3b$best[-c(1:8)]),weibull3b$best[-c(1:8)],weibull.pattern,c(rep("surv",6),rep("long",10),rep("",4)),weibull2.se)
        names(weibull.best) <- c("var","val","LC","model","se")
        weibull.best %<>% mutate(param = gsub("class\\d","",var),newvar = paste0(model," ",gsub("class\\d",paste0("class"),var),LC)) %>% arrange(model,param,LC,ascending = TRUE) %>% filter(var %ni% c('stderr','varcov 1', 'varcov 2', 'varcov 3'))
        
        weibull.comparison.data <- data.frame(weibull.best,real.params) %>% mutate(bias = val - real.params,
                                                                                   coverage = I(abs(val - real.params) < 2*se))
        
        weibull.se.estimated <- TRUE
      }
      if(is.na(weibull2.se[1])){
        weibull.best <- data.frame(names(eval(parse(text=best.weibull.model))$best[-c(1:8)]),eval(parse(text=best.weibull.model))$best[-c(1:8)],weibull.pattern,c(rep("surv",6),rep("long",10),rep("",4)),weibull.se)
        names(weibull.best) <- c("var","val","LC","model","se")
        weibull.best %<>% mutate(param = gsub("class\\d","",var),newvar = paste0(model," ",gsub("class\\d",paste0("class"),var),LC)) %>% arrange(model,param,LC,ascending = TRUE) %>% filter(var %ni% c('stderr','varcov 1', 'varcov 2', 'varcov 3'))
        
        weibull.comparison.data <- data.frame(weibull.best,real.params) %>% mutate(bias = val - real.params,
                                                                                   coverage = NA_real_)
        
        weibull.se.estimated <- FALSE
      }
    }
    
  }
  #If three-class solution not selected, then model comparison stats are set to missing.
  if(best.weibull.model != 'weibull3'){
    weibull.correct.model <- FALSE
    weibull.se <- NA_real_
    weibull.best <- NA_real_
    weibull.comparison.data <- NA_real_
    weibull.se.estimated <- NA_real_
    
  }
  
  ## Same for HARE as was done with Weibull
  if(best.hare.model == 'hare3'){
    hare.correct.model <- TRUE
    hare.se <- matrix(sqrt(diag(VarCov(eval(parse(text=best.hare.model)))))[-c(1:((hazard.k+2)*3 + 2))])
    if(!is.na(hare.se[1])){
      hare.best <- data.frame(names(eval(parse(text=best.hare.model))$best[-c(1:((hazard.k+2)*3 + 2))]),eval(parse(text=best.hare.model))$best[-c(1:((hazard.k+2)*3 + 2))],hare.pattern,c(rep("surv",6),rep("long",10),rep("",4)),hare.se)
      names(hare.best) <- c("var","val","LC","model","se")
      hare.best %<>% mutate(param = gsub("class\\d","",var),newvar = paste0(model," ",gsub("class\\d",paste0("class"),var),LC)) %>% arrange(model,param,LC,ascending = TRUE) %>% filter(var %ni% c('stderr','varcov 1', 'varcov 2', 'varcov 3'))
      
      hare.comparison.data <- data.frame(hare.best,real.params) %>% mutate(bias = val - real.params,
                                                                           coverage = I(abs(val - real.params) < 2*se))
      
      hare.se.estimated <- TRUE
    }
    if(is.na(hare.se[1])){
      hare3b <- gridsearch(rep = 50, maxiter = 25, minit = hare1,
                          Jointlcmm(y ~ visit + sex + age,
                                    random =~ 1 + visit,
                                    mixture =~ visit + sex,
                                    survival = Surv(s.inv, delta) ~ mixture(sex) + mixture(age),
                                    hazard = hazard.str,
                                    hazardnodes = knot.seq[1:knot.k],
                                    logscale = FALSE,
                                    subject = 'ID',
                                    data = data.use,
                                    partialH = FALSE,
                                    ng = 3))
      hare2.se <- matrix(sqrt(diag(VarCov(hare3b)))[-c(1:((hazard.k+2)*3 + 2))])
      if(!is.na(hare2.se[1])){
        hare.best <- data.frame(names(hare3b$best[-c(1:((hazard.k+2)*3 + 2))]),hare3b$best[-c(1:((hazard.k+2)*3 + 2))],hare.pattern,c(rep("surv",6),rep("long",10),rep("",4)),hare2.se)
        names(hare.best) <- c("var","val","LC","model","se")
        hare.best %<>% mutate(param = gsub("class\\d","",var),newvar = paste0(model," ",gsub("class\\d",paste0("class"),var),LC)) %>% arrange(model,param,LC,ascending = TRUE) %>% filter(var %ni% c('stderr','varcov 1', 'varcov 2', 'varcov 3'))
        
        hare.comparison.data <- data.frame(hare.best,real.params) %>% mutate(bias = val - real.params,
                                                                             coverage = I(abs(val - real.params) < 2*se))
        
        hare.se.estimated <- TRUE
      }
      if(is.na(hare2.se[1])){
        hare.best <- data.frame(names(eval(parse(text=best.hare.model))$best[-c(1:((hazard.k+2)*3 + 2))]),eval(parse(text=best.hare.model))$best[-c(1:((hazard.k+2)*3 + 2))],hare.pattern,c(rep("surv",6),rep("long",10),rep("",4)),hare.se)
        names(hare.best) <- c("var","val","LC","model","se")
        hare.best %<>% mutate(param = gsub("class\\d","",var),newvar = paste0(model," ",gsub("class\\d",paste0("class"),var),LC)) %>% arrange(model,param,LC,ascending = TRUE) %>% filter(var %ni% c('stderr','varcov 1', 'varcov 2', 'varcov 3'))
        
        hare.comparison.data <- data.frame(hare.best,real.params) %>% mutate(bias = val - real.params,
                                                                             coverage = NA_real_)
        
        hare.se.estimated <- FALSE
      }
    

    }
  }
  
  if(best.hare.model != 'hare3'){
   hare.correct.model <- FALSE
   hare.se <- NA_real_
   hare.best <- NA_real_
   hare.comparison.data <- NA_real_
   hare.se.estimated <- NA_real_
  }
  
  
  output_name <- paste0('Simulations_',simul_num)
  sink(paste0(your_directory,output_name,'.txt'))
  cat('###################### \n')
  cat('\n')
  cat('#PERFORMANCE \n')
  cat('\n')
  cat('###################### \n')
  cat('\n')
  cat('#Weibull Correct Model: ',weibull.correct.model,'\n')
  cat('#HARE Correct Model: ',hare.correct.model,'\n')
  cat('#Weibull SE Estimated: ',weibull.se.estimated,'\n')
  cat('#HARE SE Estimated: ',hare.se.estimated,'\n')
  cat('###################### \n')
  cat('\n')
  cat('#SUMMARY TABLE \n')
  cat('\n')
  cat('###################### \n')
  cat('\n')
  print(sum.table)
  cat('###################### \n')
  cat('\n')
  cat('#BEST WEIBULL MODEL: ',best.weibull.model,'\n')
  cat('\n')
  cat('###################### \n')
  summary(eval(parse(text=best.weibull.model)))
  cat('###################### \n')
  cat('\n')
  cat('#BEST HARE MODEL: ',best.hare.model,'\n')
  cat('\n')
  cat('###################### \n')
  cat('\n')
  summary(eval(parse(text=best.hare.model)))
  sink()
  write.csv(weibull.comparison.data,paste0(your_directory,simul_num,".csv"))
  write.csv(hare.comparison.data,paste0(your_directory,simul_num,".csv"))
  

}

#Simple trycatch function that will skip over possible replicates with errors

trycatch_jlcmm <- function(x){
  return(tryCatch(jlcmm_comparisons(x), error=function(e) NULL))
}

lapply(jlcmm.data.list,trycatch_jlcmm)




