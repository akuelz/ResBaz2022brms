
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
library(ggfortify)

# set wd
setwd("~/Desktop/Conference Materials/ResBaz2022brms/csvFiles")

########################### brms Basics (GLM) ###########################

# load data
d0 <- read.csv("dyad_cross.csv")
# Note: This is dyadic, cross-sectional data from self-identified heterosexual couples. 
  # Let's start with only one partner from each dyad to learn the basics

d1 <- subset(d0, d0$Person < 500) # removing Male partner (for now)
str(d1)
summary(d1)

# Research Question: Do individual or relational factors predict relationship conflict?

### Likelihoods ###
hist(d1$conflict) #fairly normal with some skewness
# let's compare Gaussian (Normal), Skew Normal, and Student t distributions

# Start with naive models that only include intercept and sigma
b0N <- brm(conflict ~ 1, data = d1, family = gaussian(link = "identity"), 
          chains = 4, iter = 5000, seed = 426)

b0SN <- brm(conflict ~ 1, data = d1, family = skew_normal(link = "identity"), 
            chains = 4, iter = 5000, seed = 426)

b0ST <- brm(conflict ~ 1, data = d1, family = student(link = "identity"), 
            chains = 4, iter = 5000, seed = 426)


# Note: Even very simple models can take a while to run. I always save the fitted objects so that I can read them back in
saveRDS(b0N, "simplemodelN.rds") 
saveRDS(b0SN, "simplemodelSN.rds")
saveRDS(b0ST, "simplemodelST.rds")

setwd("~/Desktop/Conference Materials/ResBaz2022brms/Saved_rdsObjects")
b0N <- readRDS("simplemodelN.rds")
b0SN <- readRDS("simplemodelSN.rds")
b0ST <- readRDS("simplemodelST.rds")

# Summary statistics of fitted models
summary(b0N)
summary(b0SN)
summary(b0ST)

# Visual check of how closely model predictions are to data
pp_check(b0N, ndraws = 50)
pp_check(b0SN, ndraws = 50)
pp_check(b0ST, ndraws = 50) # all pretty close. Skew normal appears best

### Model Comparison (cross-validation) ###
b0N <- add_criterion(b0N, c("waic", "loo", "kfold"))
b0SN <- add_criterion(b0SN, c("waic", "loo", "kfold"))
b0ST <- add_criterion(b0ST, c("waic", "loo", "kfold"))

loo_compare(b0N, b0SN, b0ST, criterion="waic") 
loo_compare(b0N, b0SN, b0ST, criterion="loo") 
loo_compare(b0N, b0SN, b0ST, criterion="kfold") 
# elpd_diff is the "expected log predictive density" which is the difference between the models in terms of expected predictive accuracy in new data
# first line compares the model with the highest elpd to itself, so the values are always zero
# next line(s) compare the best model to the others
# consider the diff relative to its SE to decide what the evidence is

### Model Estimation ###
b1 <- brm(conflict ~ relstressA + esteem, data = d1, family = skew_normal(link = "identity"), 
          chains = 4, iter = 5000, seed = 426 )
b1 <- add_criterion(b1, criterion = c("waic", "loo", "kfold"))
saveRDS(b1, "b1_GLM.rds") 
b1 <- readRDS("b1_GLM.rds")

summary(b1) # summary statistics of fitted model
plot(b1) # visualize full posterior distribution and chains
conditional_effects(b1)
bayes_R2(b1)

### Posterior Predictive Checks ###

# Visual check of how closely model predictions are to data
pp_check(b1, ndraws = 50)

# Check Residuals
pp_check(b1, type = "error_hist", ndraws = 20)

# Homoscedasticity across values of Y
pp_check(b1, type="error_scatter_avg", ndraws=50) 

# Homoscedasticity across values of X
pp_check(b1, type = "error_scatter_avg_vs_x", x="relstressA", ndraws = 20)
pp_check(b1, type = "error_scatter_avg_vs_x", x="esteem", ndraws = 20)

# What other pp checks are there? 
pp_check(b1N, type = "xyz")

### ROPEs ###
rope(b1) # computes the proportion (in percentage) of the HDI (default to the 89% HDI) of a posterior distribution that lies within a region of practical equivalence.
rope(b1, range= rope_range(b1), ci=0.95)  # sets the range to 0 +/- .1 * sd(y) and the HDI to 95%
plot(rope(b1, range = rope_range(b1), ci = 0.95)) # plot the rope
equivalence_test(b1)

### Priors ###

# First find out what default priors are. 
prior_summary(b1)
# The student-t priors have 3 parameters: degrees of freedom, mu, sigma. 
  # The Student’s t distribution that emerges for the y-intercept is going to have three degrees of freedom, a mu that equals the median of the dependent variable, and that sigma it selects. 
  # The residual standard deviation is going to have a Student’s t with three degrees of freedom, a mu of zero, and that sigma it selects.

median(d1$conflict, na.rm = T)
mad(d1$conflict, na.rm = T)
# If the median absolute deviation of the dependent variable is less than 2.5, the sigma it chooses is going to be 2.5. 
# If the median absolute deviation of the dependent variable is greater than 2.5, it will round that to a single decimal point and use that as the sigma. 
 
## We can look at these visually
# Intercept: plot a t that is centered around mu=2 and goes from +- 3 SDs and compare it to an uninformative uniform distribution (the default is +1 INF, so we just use big numbers for the range)
p <- ggdistribution(dstudent_t, seq(-6, 10, .05), df=3, mu=2, sigma=2.5, colour='blue')
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

b2 <- brm(conflict ~ relstressA + esteem, data=d1, prior = ourPrior, family=skew_normal(link = "identity"), 
          chains = 4, iter = 5000, seed = 426)
b2 <- add_criterion(b2, c("waic", "loo", "kfold"))
saveRDS(b2, "b2_GLM.rds")
b2 <- readRDS("b2_GLM.rds")
prior_summary(b2)
summary(b2) # estimate is shrunk a bit, but we still reach the same conclusions
bayes_R2(b2)

loo_compare(b1, b2, criterion="loo") 
loo_compare(b1, b2, criterion = "kfold")

########################### Complex Models ###########################

# Remember that this was dyadic data? 

### Cross-Sectional Dyadic Model ###
hist(d0$conflict) # let's go with skew normal likelihood

# Random-Intercept model with no predictors 
b3 <- brm(conflict ~ 1 + (1|dyad), data=d0, family=skew_normal(link = "identity"), 
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
icc1 # 38% of variance due to between dyad differences; 62% due to differences within dyad

## actor-partner with random intercept model
b4 <- brm(conflict ~ relstressA + p_relstressA + (1|dyad), data=d0, control= list(adapt_delta=.99), 
          family=skew_normal(link = "identity"), chains = 4, iter = 5000, seed = 426)
# note that this is a relatively simple model, but there were divergent transitions. so, we increased adapt_delta
saveRDS(b4, "crossSectional_b4.rds")
b4 <- readRDS("crossSectional_b4.rds")
summary(b4)
conditional_effects(b4)
pp_check(b4, type="error_hist", ndraws=20)
pp_check(b4, type = "error_scatter_avg_vs_x", x="relstressA", ndraws = 50)
pp_check(b4, type = "error_scatter_avg_vs_x", x="p_relstressA", ndraws = 50)


bayes_R2(b3)
bayes_R2(b4)


### Repeated Measures Dyadic Model ###
setwd("~/Desktop/Conference Materials/ResBaz2022brms/csvFiles")
dRep <- read.csv("dyad_repeatWide.csv")
str(dRep)
summary(dRep)

# time nested within person within dyad
# Bayesian multivariate models contain separate parts for each outcome and need the data in wide format. 
  # We will just use default likelihood and priors to keep the focus on the models.

hist(dRep$f_hb)
hist(dRep$m_hb)

## Two-intercept model (with no predictors), but full error structure

bf_f <- bf(f_hb~ 1 + (1|p|DYAD), autocor= ~ar(time= DAY, gr=DYAD)) 
# |p| indicates that all varying effects of DYAD should be modeled as correlated. 
  # The indicator p is arbitrary and can be replaced by any other symbol that comes to mind 
bf_m <- bf(m_hb~ 1 + (1|p|DYAD), autocor= ~ar(time= DAY, gr=DYAD))
b5 <- brm(bf_f + bf_m + set_rescor(TRUE), data=dRep, chains=2, cores=10, seed = 426)
saveRDS(b5, "repeated_b5.rds")

setwd("~/Desktop/Conference Materials/ResBaz2022brms/Saved_rdsObjects")
b5 <- readRDS("repeated_b5.rds")
summary(b5)

## Growth curves, with random intercepts but fixed slopes. Also add fixed effects of actor and partner closeness for the sake of an example.
bf_f <- bf(f_hb~ DAY + f_CLOSE + m_CLOSE + (1|p|DYAD), autocor= ~ar(time= DAY, gr=DYAD))
bf_m <- bf(m_hb~ DAY + f_CLOSE + m_CLOSE + (1|p|DYAD), autocor= ~ar(time= DAY, gr=DYAD))
b6 <- brm(bf_f + bf_m + set_rescor(TRUE), data=dRep, chains=2, cores=10, seed = 426)
saveRDS(b6, "repeated_b6.rds")
b6 <- readRDS("repeated_b6.rds")
summary(b6)

## And just to show the versatility, here we have completely different models for the women and men, 
  # including different fixed effects, autocorrelation for the men but not women and a 
  # skew-normal likelihood for the women but the default Gaussian for the men. 
  # Note: We can no longer have the residual correlation between men and women if we use different likelihoods.
bf_f <- bf(f_hb~  f_CLOSE + (1|p|DYAD), family="skew_normal")
bf_m <- bf(m_hb~  m_ARGUE + (1|p|DYAD), autocor= ~ar(time= DAY, gr=DYAD))
b7 <- brm(bf_f + bf_m, data=dRep, chains=2, cores=10, seed = 426)
saveRDS(b7, "repeated_b7.rds")
b7 <- readRDS("repeated_b7.rds")
summary(b7)

