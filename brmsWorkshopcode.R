
############################
### Introduction to brms ###
############################

# install necessary packages
#install.packages(brms)
#install.packages("bayestestR")

# load packages
library(brms)
library(bayestestR)
library(bayesplot)

# set wd
#setwd("~/Desktop/Conference Materials/ResBaz")

########################### brms Basics (GLM) ###########################

# load data
d1 <- read.csv("regressionExample.csv", stringsAsFactors = T)

str(d1)
summary(d1)

### Model Estimation ###
b1N <- brm(conflict ~ relstressA, data = d1, family = gaussian(link = "identity"), 
          chains = 4, iter = 5000, seed = 426)

# models can take a while to run. I always save the fitted object so that I can read it back in
saveRDS(b1N, "simplemodelN.rds") 
b1N <- readRDS("simplemodelN.rds")

summary(b1N) # summary statistics of fitted model
plot(b1N) # visualize full posterior distribution and chains
conditional_effects(b1N)

### Likelihoods ###
hist(d1$conflict) # fairly normal, but somewhat skewed. We will compare a normal likelihood (Gaussian) to skewed-normal and student_t likelihoods

b1ST <- brm(conflict ~ relstressA, data = d1, family = student(link = "identity"), 
            chains = 4, iter = 5000, seed = 426)
b1SN <- brm(conflict ~ relstressA, data = d1, family = skew_normal(link = "identity"), 
            chains = 4, iter = 5000, seed = 426)

saveRDS(b1ST, "simplemodelST.rds")
saveRDS(b1SN, "simplemodelSN.rds")

b1ST <- readRDS("simplemodelST.rds")
b1SN <- readRDS("simplemodelSN.rds")

### Posterior Predictive Checks ###

# Visual check of how closely model predictions are to data
pp_check(b1N, ndraws = 50)
pp_check(b1ST, ndraws = 50)
pp_check(b1SN, ndraws = 50)

# Check Residuals
pp_check(b1N, type = "error_hist", ndraws = 20)
pp_check(b1ST, type = "error_hist", ndraws = 20)
pp_check(b1SN, type = "error_hist", ndraws = 20)

# Homoscedasticity across values of Y
pp_check(b1N, type="error_scatter_avg", ndraws=50) 
pp_check(b1ST, type="error_scatter_avg", ndraws=50)
pp_check(b1SN, type="error_scatter_avg", ndraws=50)

# Homoscedasticity across values of X
pp_check(b1N, type = "error_scatter_avg_vs_x", x="relstressA", ndraws = 20)
pp_check(b1ST, type = "error_scatter_avg_vs_x", x="relstressA", ndraws = 20)
pp_check(b1SN, type = "error_scatter_avg_vs_x", x="relstressA", ndraws = 20)

# What other pp checks are there? 
pp_check(b1N, type = "xyz")

### Model Comparison (cross-validation) ###
b1N <- add_criterion(b1N, c("waic", "loo", "kfold"))
b1ST <- add_criterion(b1ST, c("waic", "loo", "kfold"))
b1SN <- add_criterion(b1SN, c("waic", "loo", "kfold"))

loo_compare(b1N, b1ST, b1SN, criterion="waic") 
loo_compare(b1N, b1ST, b1SN, criterion="loo") 
loo_compare(b1N, b1ST, b1SN, criterion="kfold") 
# elpd_diff is the "expected log predictive density" which is the difference between the models in terms of expected predictive accuracy in new data
# first line compares the model with the highest elpd to itself, so the values are always zero
# next line(s) compare the best model to the others
# consider the diff relative to its SE to decide what the evidence is

### ROPEs ###
rope(b1SN) # computes the proportion (in percentage) of the HDI (default to the 89% HDI) of a posterior distribution that lies within a region of practical equivalence.
rope(b1SN, range= rope_range(b1SN), ci=0.95)  # sets the range to 0 +/- .1 * sd(y) and the HDI to 95%
plot(rope(b1SN, range = rope_range(b1SN), ci = 0.95)) # plot the rope
equivalence_test(b1SN)

### Priors ###

# First find out what default priors are. 
prior_summary(b1SN)
# The student-t priors have 3 parameters: degrees of freedom, mu, sigma. 
  # The Student’s t distribution that emerges for the y-intercept is going to have three degrees of freedom, a mu that equals the median of the dependent variable, and that sigma it selects. 
  # The residual standard deviation is going to have a Student’s t with three degrees of freedom, a mu of zero, and that sigma it selects.

median(d1$conflict)
mad(d1$conflict)
# If the median absolute deviation of the dependent variable is less than 2.5, the sigma it chooses is going to be 2.5. 
# If the median absolute deviation of the dependent variable is greater than 2.5, it will round that to a single decimal point and use that as the sigma. 
 
## We can look at these visually
# Intercept: plot a t that is centered around mu=1.8 and goes from +- 3 SDs and compare it to an uninformative uniform distribution (the default is +1 INF, so we just use big numbers for the range)
p <- ggdistribution(dstudent_t, seq(-6, 10, .05), df=3, mu=1.8, sigma=2.5, colour='blue')
p <- ggdistribution(dunif, seq(-6, 10, .05), min=-1000, max=1000, colour="red", p=p)
p

# Sigma: plot a t that is centered around mu=0 and goes from +- 3 SDs and compare it to an uninformative uniform distribution (the default is +1 INF, so we just use big numbers for the range)
p <- ggdistribution(dstudent_t, seq(-8, 8, .05), df=3, mu=0, sigma=2.5, colour='blue')
p <- ggdistribution(dunif, seq(-8, 8, .05), min=-1000, max=1000, colour="red", p=p)
p

# Let's try two different priors: one weakly informative prior and one strong prior
# Let's say there was good prior research suggesting that the beta for relstressA should be about 0.5 (e.g., smaller than we found so far)

### figure out what distribution would describe an effect centered at 0.5, but ranging from about 0 to 1

x <- seq(-1, 2, .1) # create sequence for x-axis 
plot(x, dnorm(x, mean= 0.5, sd=.25), type="l", ylim=c(0,2)) # plot the normal distribution we specified

## to get a weaker prior (e.g., more spread out) we can increase the sd
plot(x, dnorm(x, mean= 0.5, sd=.5), type="l", ylim=c(0,2)) 

# we'll go with the fairly strong prior
ourPrior <- set_prior("normal(0.5,.25)", class="b", coef = "relstressA")

weakPrior <- set_prior("uniform(-10,10)", class = "b",  lb = -10, ub = 10)

b2 <- brm(conflict ~ relstressA , data=d1, prior = ourPrior, family=skew_normal(link = "identity"), 
          chains = 4, iter = 5000, seed = 426)
b2 <- add_criterion(b2, c("waic", "loo", "kfold"))
saveRDS(b2, "simpleModelPrior.rds")
b2 <- readRDS("simpleModelPrior.rds")
prior_summary(b2)
summary(b2) # estimate is shrunk a bit

loo_compare(b1SN, b2, criterion="loo") 
loo_compare(b1SN, b2, criterion = "kfold")

bayes_R2(b1SN)
bayes_R2(b2)



m1SN <- brm(posRel ~ pospexpress, data = d1, family = skew_normal, chains = 4, iter = 5000, seed = 426)
saveRDS(m1SN, "simplemodelSN.rds")
m1SN <- readRDS("simplemodelSN.rds")
summary(m1SN)

pp_check(m1SN, ndraws = 50)
pp_check(m1N, ndraws = 50)


########################### Complex Models ###########################

### Cross-Sectional Dyadic Model ###
dCross <- read.csv("dyad_cross.csv")

hist(dCross2$stressTot) # let's go with skew normal likelihood

# Random-Intercept model with no predictors 
b3 <- brm(stressTot ~ 1 + (1|dyad), data=dCross, family=skew_normal(link = "identity"), 
            chains = 4, iter = 5000, seed = 426)
saveRDS(b3, "crossSectional_b3.rds")
b3 <- readRDS("crossSectional_b3.rds")
summary(b3)
pp_check(b3, ndraws=50)

# ICC
iSd <- as.numeric(summary(b3)$random$dyad[1,1])
iVar <- iSd^2
residSd <- as.numeric(summary(b3)$spec_pars[[1,1]])
residVar <- residSd^2
icc1 <- iVar/(iVar+residVar) 
icc1 # 27% of variance due to between dyad differences; 73% due to differences within dyad

## actor-partner with random intercept model
b4 <- brm(stressTot ~ conflict + p_conflict + (1|dyad), data=dCross, control= list(adapt_delta=.99), 
          family=skew_normal(link = "identity"), chains = 4, iter = 5000, seed = 426)
# note that this is a relatively simple model, but there were divergent transitions. so, we increased adapt_delta
saveRDS(b4, "crossSectional_b4.rds")
b4 <- readRDS("crossSectional_b4.rds")
summary(b4)

bayes_R2(b3)
bayes_R2(b4)


### Repeated Measures Dyadic Model ###
dRep <- read.csv("dyad_repeatWide.csv")

# time nested within person within dyad
# Bayesian multivariate models contain separate parts for each outcome and need the data in wide format. 
  # We will just use default likelihood and priors to keep the focus on the models.

hist(dRep$f_hb)
hist(dRep$m_hb)

## Two-intercept model (with no predictors), but full error structure; note we now get separate autocorrelation estimates for men and women, as well as everything else

bf_f <- bf(f_hb~ 1 + (1|p|DYAD), autocor= ~ar(time= DAY, gr=DYAD)) # |p| indicates that all varying effects of DYAD should be modeled as correlated. The indicator p is arbitrary and can be replaced by other symbol that comes to mind 
bf_m <- bf(m_hb~ 1 + (1|p|DYAD), autocor= ~ar(time= DAY, gr=DYAD))
b1_dyad <- brm(bf_f + bf_m + set_rescor(TRUE), data=data, chains=2, cores=10)
saveRDS(b1_dyad, "repeated2_b1_dyad.rds")
b1_dyad <- readRDS("repeated2_b1_dyad.rds")
summary(b1_dyad)


b1SN <- brm(hb ~ 1 + (1|Person), data=d2, family=skew_normal(link = "identity"), 
            chains = 4, iter = 5000, seed = 426)
b1SN <- add_criterion(b1SN, c("waic", "loo", "kfold"))
saveRDS(b1SN, "uncondModel.rds")

hist(d2$MeanSC)
hist(d2$IBI)

# Unconditional means model
b1SN <- brm(MeanSC ~ 1 + (1|Person), data=d2, family=skew_normal(link = "identity"), 
            chains = 4, iter = 5000, seed = 426)
b1SN <- add_criterion(b1SN, c("waic", "loo", "kfold"))
saveRDS(b1SN, "uncondModel.rds")



# Time in variable in person
d2 <- subset(d2, time < 25)

# Growth curves with fixed slopes but  otherwise full error structure
bf_sc <- bf(MeanSC~ time + (1|p|Person), autocor= ~ar(time= time, gr=Person))
bf_ibi <- bf(IBI~ time + (1|p|Person), autocor= ~ar(time= time, gr=Person))
b1_bivar <- brm(bf_sc + bf_ibi + set_rescor(TRUE), data=d2, chains=2, cores=10)
saveRDS(b1_bivar, "repeated2_b1_bivar.rds")
b1_bivar <- readRDS("repeated2_b1_bivar.rds")
summary(b1_bivar)


d2 <- read.csv("diaryClean.csv", stringsAsFactors = T)
# diary data (2 reports a day for 10 days)

# Unconditional means model
b1SN <- brm(posRel ~ 1 + (1|Person), data=d2, family=skew_normal(link = "identity"), 
            chains = 4, iter = 5000, seed = 426)
b1SN <- add_criterion(b1SN, c("waic", "loo", "kfold"))
saveRDS(b1SN, "uncondModel.rds")


b1SN <- readRDS("uncondModel.rds")
summary(b1SN)

# Check ICC to get sense of how variance is distributed between- and within- people
iSd <- as.numeric(summary(b1SN)$random$Person[1,1])
iVar <- iSd^2
residSd <- as.numeric(summary(b1SN)$spec_pars[[1,1]])
residVar <- residSd^2
ICCb <- iVar/(iVar+residVar)
ICCb # ~ 17% due to between person differences, remaining 83% is due to within-person differences over time

b2 <- brm(posRel ~ 1 + time + pospexpress + t1ambiv + (1 + time + pospexpress|Person), data=d2, family=skew_normal(link="identity"), 
          control=list(adapt_delta = .95), chains = 4, iter = 5000, seed = 426)
saveRDS(b2, "fullModel.rds")
summary(b2)
