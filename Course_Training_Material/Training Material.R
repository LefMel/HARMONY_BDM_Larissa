# HARMONY CA18208 Training School 

# Introduction to Bayesian Disease Measurement (BDM) for health scientists


# 25 - 27 October 2022
# Larissa, Greece


# Day 1 Material (Polychronis Kostoulas, Eleftherios Meletis)
#################################

# Bayes theorem - Predictive Values / Exercise

# A rapid test has been developed to detect if a person is infected with the new SARS-CoV-2 virus. This test is fairly reliable:

#  * 95% of all infected individuals are detected and, 
#  * 95% of all healthy individuals are identified as such. 
#  * Also, it has been documented that at most one passenger out of 100 aboard on an airplane is infected.

#  A1. What is the Sensitivity and Specificity of the test and the prevalence of the population?
#  A2. Estimate the probability of of a person being infected given that he/she tested positive?
#  A3. Estimate the probability of a person being healthy, given that he/she tested negative?

# Alterinative case - What if:

#  * 60% of all infected individuals are detected and,
#  * 99% of all healthy individuals are identified as such, the Sp of the test is 100%
#  * The prevalence can be 100% or 0%.

#  B1. Estimate now the predictive values of the test.

#  C1. What are the predictive values when we have a "perfect" test?




# Solution
# A1. Sensitivity = 95%, Specificity = 95%, prevalence = 1 / 100 = 1%
# A2. The metric we want to calculate is the Positive Predictive value (PPV), using Bayes theorem
# A3. The metric we want to calculate is the Negative Predictive value (NPV), using Bayes theorem

# Let's estimate the PPV and NPV with the help of a R-package (BioProbability)
# R - package to estimate PVs
install.packages("BioProbability")
library(BioProbability)
p<-seq(0.01,0.99, by=0.01) # Range of values (0-100) for the prevalence
predictive.value(p,Spe=0.95,Sen=0.95,plot.it=TRUE)

# prev = 0.01, Spe=0.95, Sen = 0.95
predictive.value(0.01,Spe=0.95,Sen=0.95)
#         Prevalence + predictive value - predictive value
#[1,]       0.01          0.1610169          0.9994687

# B1
# prev = p, Spe=0.99, Sen = 0.6
predictive.value(p,Spe=0.99,Sen=0.6, plot.it=TRUE)

# C1  
# prev = p, Spe=1, Sen=1
predictive.value(p,Spe=1,Sen=1, plot.it=TRUE)


# Topics

### During this school the following topics will be covered:

#  1. Running basic models in JAGS
#  2. Apparent & true prevalence estimation
#  3. How to choose a prior distribution?
#  4. Diagnostic test evaluation with Hui-Walter models

# Main questions

### What is JAGS and what is MCMC?

### How are the apparent and true prevalence defined and why are they different?

### What is a prior distribution?

### How can we estimate the Sensitivity and Specificity of a diagnostic test?


# Probability distributions

# * Likelihood theory is at the heart of most inferential statistics

# A likelihood is the probability of observing our data given the distribution that we use to describe the data generating process.  

# Example

# * What is the likelihood (i.e. probability) of getting 5 heads from 10 tosses of a fair coin?  

# We assume:

# * The probability distribution describing a number of independent coin tosses is called the Binomial distribution
# * In this case, we would use the parameters:
#	+ Number of coin tosses = 10
#	+ Probability of a head = 'fair' = 0.5

# Example - Rstudio
tosses <- 10
probability <- 0.5
heads <- 5
likelihood_1 <- choose(tosses, heads) * probability^heads *
  (1-probability)^(tosses-heads)
likelihood_1
# [1] 0.2460938


#But R makes our life easier by implementing this using a function called dbinom:
?dbinom()
likelihood_2 <- dbinom(heads, tosses, probability)
likelihood_2
# [1] 0.2460938


# Maximising a likelihood

# In the previous example we assumed that we knew the probability of getting a head because the coin was fair (i.e. probability of head = 50%), but typically we would want to estimate this parameter based on the data.  

# One way to do this is via Maximum Likelihood.  

# Let's say that we have observed 7 test positive results from 10 individuals and we want to estimate the prevalence by maximum likelihood.  

# We could do that by defining a function that takes our parameter as an argument, then calculates the likelihood of the data based on this parameter:
likelihood_fun <- function(prevalence) dbinom(7, 10, prevalence)

#We can now ask the function what the likelihood is for any parameter value that we choose, for example:
likelihood_fun(0.8)
# [1] 0.2013266 
likelihood_fun(0.5)
# [1] 0.1171875

# So the data are more likely if the prevalence parameter is 0.8 than if it is 0.5.  

# We could keep doing this for lots of different parameter values until we find the highest likelihood, but it is faster and more robust to use an R function called optimise to do this for us:
optimise(likelihood_fun, interval=c(0, 1), maximum=TRUE)
# $maximum
# [1] 0.6999843

# $objective
# [1] 0.2668279

# This tells us that the maximum likelihood for this data is 0.267, which corresponds to a parameter value of around 0.7 (or a prevalence of 70%).  This is the maximum likelihood.  


# Profiling a likelihood

# The parameters corresponding to the maximum likelihood give the highest probability of observing the data given the parameters, but there are other parameter values under which we could observe the data with almost as high a probability.  

# It is useful to look at the range of parameter values that are consistent with the data, which is why R reports standard errors (and/or confidence intervals) when you run a model. 

# But we can also look at the full distribution of the likelihood of the data over a range of parameter values using our function above.

# Example
parameters <- seq(0, 1, length.out=101)
likelihoods <- numeric(length(parameters))
for(i in 1:length(parameters)){
  likelihoods[i] <- likelihood_fun(parameters[i])
}
plot(parameters, likelihoods, type='l')
abline(h=0.267, lty='dashed', col='red')
abline(v=0.7, lty='dashed', col='red')

#The red dashed lines show the maximum likelihood (y axis) with corresponding parameter value (x axis), and the solid line is the likelihood of the data given the parameter value on the x axis.  You can see that parameter values near 0.7 have a likelihood that is almost as high as the maximum.

# Bayesian Statistics

#In this session we'll see how we can estimate a probability of interest  but in a Bayesian framework, i.e. using Bayes theorem.

# Bayes' theorem 

#  P(A|B) = P(B|A)*P(A)/P(B)


# Components

# * P(A|B): Prob of event A occurring given that B is true - Posterior probability
# * P(B|A): Prob of event B occurring given that A is true - Likelihood ~ function of A
# * P(A): Prob of event A occurring - Prior probability
# * P(B): Prob of event B occurring - Marginal probability ~ sum over all possible values of A

# What we usually see/use


# theta: parameter of interest | y: observed data}
# P(theta|y) = P(y|\theta) * P(theta)/P(y) 
# Where:

#  * P(theta): Prior probability of parameter(s) of interest;
#  * P(y|theta): Likelihood of the data given the parameters value(s) 
#  * P(theta|y): Posterior probability of parameter(s) of interest given the data and the prior

# Bayesian Inference - Summary & Example

# To estimate the posterior distribution P(theta$|y) we need to:

#  Specify the Prior distribution: P(theta)
#  Define the Likelihood of the data: P(y|theta) 

# Example: Bayesian apparent prevalence (ap) estimation

# y out of n individuals test positive. Estimate the apparent prevalence.


# Parameter of interest: ap - [0,1]

# Data: n tested, y positive

# * Prior distribution for ap: ap ~ Beta(a,b)
# * Likelihood: y ~ Binomial(n,ap) 

# Let's write our first JAGS model

ap_model <- 
  'model {
  
  # Define likelihood distribution of the data
  # JAGS Binomial distribution Arguments: p, n 
  
  y ~ dbin(ap,n)
  
  # Specify prior distribution for par of interest 
  # Uniform (non-informative) prior distribution 
  ap ~ dbeta(1, 1)

  #data# n, y
  #monitor# ap
  #inits# ap
  }
  '


# Let's run our first JAGS model

# Call JAGS
library(runjags)

# Provide Data 
n = 40
y = 12

# Initial values for par of interest
ap <- list(chain1=0.05, chain2=0.95)


# Run the model
results <- run.jags(ap_model, n.chains=2, 
                    burnin=5000, sample=10000)

# View results

# Plot results
plot(results)
# Print results
summary(results)

#################################

# Day 2 Material (Konstantinos Pateras, Julio Alvarez, Matthew Denwood)
#################################

# How to choose and generate prior distributions?

# Kostas Pateras pdfs on PriorGen and tPriors available in GitHub - Course_training_material folder

require(PriorGen)

# Exercise 1: Check if the codes below create the prior parameters for slides 1 to 3.

# Example 1
findbeta(themean = .9, percentile = .95, lower.v = F, percentile.value = .80)
# Example 2
findbeta(themean = .99, percentile = .95, lower.v = F, percentile.value = .90)
# Example 3
findbeta(themean = .01, percentile = .70, lower.v = T, percentile.value = .05)


# Exercise 2: What if I wanted to create a beta prior for the specificity and I knew that 
# 1. the mean specificity lies between (0.4-0.6) and 
# 2. we are 95% certain that it is lower than 0.8
findbeta(themean=0.5, percentile = 0.95, lower.v = T, percentile.value = 0.8)

# Exercise 3: Plotting of priors based on hyperparameters of Beta distribution.
x=seq(0,1,0.001)
plot(x,dbeta(x,shape1 = 27.79, 3.09),type="l",lwd=3,ylab="Density (Sensitivity)")
lines(x,dbeta(x,shape1 = 7.59, .08),type="l",lwd=3,ylab="Density (Sensitivity)", col="red")

#  Plot the density for Example 3!
x=seq(0,1,0.001)
plot(x,dbeta(x,shape1 = 3.26, shape2 = 3.26),type="l",lwd=3,ylab="Density (Sensitivity)")


# Matthew Denwood presentation

  ## An (extremely brief) introduction to Bayesian Markov chain Monte Carlo ##
      ## Accessible in GitHub and here (Ctl + click): https://drive.google.com/file/d/1na5KVYQ17vnMR3_Gr0bTONTzODdooIbr/view?usp=share_link


# Session Exercise

## Posterior function:

# Define a function that calculates a posterior:
log_posterior_fun <- function(parameter){
  # ... Some R code e.g. ...
  ll <- dbinom(data$Pos, data$N, parameter, log=TRUE)
  lp <- dbeta(parameter, 1, 1, log=TRUE)
  return(ll + lp)
}

# We can use this with the following data:
data <- list(Pos = 1210, N = 4072)

# For example:
log_posterior_fun(qlogis(0.25))

## Basic MCMC algorithm

# Set up the parameter and posterior vector:
iters <- 1000
parameter <- numeric(iters)
log_post <- numeric(iters)

# Initial values:
parameter[1] <- 0.25
log_post[1] <- log_posterior_fun(parameter[1])

# Pick a value for sigma (this can be made bigger or smaller):
sigma <- 0.01

# Metropolis algorithm loop:
for(i in 2:iters){
  
  # Sample a new parameter value:
  new_par <- rnorm(1, mean=parameter[i-1], sd=sigma)
  
  if(new_par < 0 || new_par > 1){
    # If the new parameter is invalid then reject it:
    accept <- 0
  }else{
    new_lpost <- log_posterior_fun(new_par)
    if(new_lpost > log_post[i-1]){
      # If this is an improvement always accept:
      accept <- 1
    }else{
      # Otherwise do an accept/reject step:
      probability_ratio <- exp(new_lpost - log_post[i-1])
      # This is the same as:
      # probability_ratio <- exp(new_lpost) / exp(log_post[i-1])
      accept <- rbinom(1, 1, probability_ratio)
    }
  }
  
  # Save the new value if accepted, otherwise copy the old value:
  if(accept==1){
    parameter[i] <- new_par
    log_post[i] <- new_lpost
  }else{
    parameter[i] <- parameter[i-1]
    log_post[i] <- log_post[i-1]
  }		
}

# Trace plot:
plot(parameter, type="l")

# Removing burnin:
plot(51:iters, parameter[-(1:50)], type="l", xlim=c(0,iters))

# Histogram:
hist(parameter[-(1:50)])

# Effective sample size:
library("coda")
effectiveSize(parameter[-(1:50)])
# Note: this is not enough, we need more iterations!

# Median estimate and 95% CI:
median(parameter[-(1:50)])
HPDinterval(as.mcmc(parameter[-(1:50)]))

## Equivalent in JAGS
bugs_model <- "
	model{
		# Likelihood:
		Pos ~ dbinom(ap, N)
		# Uniform (non-informative) prior for apparent prevalence
		ap ~ dbeta(1,1)
		
		#data# Pos, N
		#inits# ap
		#monitor# ap
	}
  "

# Data:
Pos <- 1210
N <- 4072

# Initial values:
ap <- list(chain1=0.1, chain2=0.9)

library('runjags')
results <- run.jags(bugs_model, n.chains=2, burnin=5000, sample=10000)

plot(results)
results

# Session - How the priors affect the posteriors? -

#initial scenario with small sample size: results from the experiment
Positives <- 35
N <- 121

ap_prev = "model{

    # Likelihood part:
    Positives ~ dbinom(prevalence, N)

    # Prior part:
    prevalence ~ dbeta(1, 1)

    # Hooks for automatic integration with R:
    #data# Positives, N
    #monitor# prevalence
    #inits# prevalence
  }
  "

# Run the following scenarios and compare the results.

# Scenario 1: non-informative (weakly informative) priors
prevalence <- list(chain1=0.05, chain2=0.95)

results1 <- run.jags(ap_prev, n.chains=2, burnin=5000, sample=10000)
plot(results1)
results1


# Scenario 2: highly informative priors: 95% sure that prevalence is below 3%, most likely value is around 0.5%
findbeta(themode = 0.05, percentile = 0.95, lower.v = F, percentile.value = 0.03) # Be(9.64, 165.2)

# We need to change the prior in our model.
ap_prev_2 = "model{

    # Likelihood part:
    Positives ~ dbinom(prevalence, N)

    # Prior part:
    prevalence ~ dbeta(9.64, 165.2)

    # Hooks for automatic integration with R:
    #data# Positives, N
    #monitor# prevalence
    #inits# prevalence
  }
  "

prevalence <- list(chain1=0.05, chain2=0.95)

results2 <- run.jags(ap_prev_2, n.chains=2, burnin=5000, sample=10000)
plot(results2)
results2


# Scenario 3: somewhat informative priors: 95% sure that prevalence is below 10%, most likely value is around 3%
findbeta(themode = 0.03, percentile = 0.95, lower.v = T, percentile.value = 0.1) # Be(2.63, 53.58)

# We need to change the prior in our model.
ap_prev_3 = "model{

    # Likelihood part:
    Positives ~ dbinom(prevalence, N)

    # Prior part:
    prevalence ~ dbeta(2.63, 53.58)

    # Hooks for automatic integration with R:
    #data# Positives, N
    #monitor# prevalence
    #inits# prevalence
  }
  "

prevalence <- list(chain1=0.05, chain2=0.95)

results3 <- run.jags(ap_prev_3, n.chains=2, burnin=5000, sample=10000)
plot(results3)
results3

# Try the 3 scenarios above but with a large sample size. 
# We do not need to change the model, only the data (Positive and N) in the R environment
# Lets try:
Positives <- 3500
N <- 12100

# Scenario 1 Larger size
results_S_1 = run.jags(ap_prev, n.chains=2, burnin=5000, sample=10000)
plot(results_S_1)
results_S_1

# Scenario 2 Larger size
results_S_2 = run.jags(ap_prev_2, n.chains=2, burnin=5000, sample=10000)
plot(results_S_2)
results_S_2

  
# Scenario 3 Larger size
results_S_3 = run.jags(ap_prev_3, n.chains=2, burnin=5000, sample=10000)
plot(results_S_3)
results_S_3

Results_Table = rbind(results1$summaries[,1:3], results2$summaries[,1:3], results3$summaries[,1:3],  results_S_1$summaries[,1:3], results_S_2$summaries[,1:3], results_S_3$summaries[,1:3])
row.names(Results_Table) = c("Small Size + Weak", "Small Size + High", "Small Size + Intermediate", "Large Size + Weak", "Large Size + High", "Large Size + Intermediate")
Results_Table

# Example: Bayesian true prevalence (tp) estimation

# Assuming the absence of a perfect test we do not know how many individuals are truly positive/negative.

# Instead we know that n individuals are tested with an imperfect test and y have a positive result.


# Apparent/True prevalence: ap/tp 
# Sensitivity: Se 
# Specificity: Sp

# ap = P(T+) = P(T + & D+) + P(T+ &  D-) = P(D+) * P(T+|D+) + P(D-) * P(T+|D-) =>
# ap = tp * Se + (1 - tp) * (1 - Sp)

## Create a JAGS model for true prevalence estimation

# Parameters of interest - tp, Se, Sp 

# Prior distributions

# tp ~ dbeta(1,1) # Uniform (non-informative) prior distribution 
# Se ~ dbeta(25.4, 3.4) # 0.85 (0.7 - 0.95)
# Sp ~ dbeta(95, 5) # 0.77 (0.49 - 0.96)

# Data: n tested, y positive

# Likelihood: y ~ Binomial(n,ap), 
# ap = tp * Se + (1 - tp) * (1 - Sp)

## Write JAGS model
tp_model <- 
  'model {
  y ~ dbin(ap,n)
  ap <- tp*Se + (1-tp)*(1-Sp)
  
  # Uniform (non-informative) prior distribution 
  tp ~ dbeta(1,1)
  # Informative priors for Se and Sp
  Se ~ dbeta(25.4, 3.4)
  Sp ~ dbeta(95, 5)
  
  #data# n, y
  #monitor# tp, Se, Sp
  #inits# tp, Se, Sp
  }
  '


## Let's run our JAGS model
# Call JAGS
library(runjags)
# Provide Data 
n = 4072
y = 1210
# Initial values for pars of interest
tp <- list(chain1=0.05, chain2=0.95)
Se <- list(chain1=0.05, chain2=0.95)
Sp <- list(chain1=0.05, chain2=0.95)

# Run the model
results <- run.jags(tp_model, n.chains=2, 
                    burnin=5000, sample=10000)
## View results
# Plot results
plot(results) # Click backwards to view all plots
# Print results
summary(results)



#################################

# Day 3 Material (Eleftherios Meletis, Julio Alvarez)
#################################

# Sensitivity - Specificity estimation with and without a gold standard

## Hui-Walter paradigm/model (1980) 

# Link to Publication: https://doi.org/10.2307/2530508 
# Estimating the Error Rates of Diagnostic Tests S. L. Hui and S. D. Walter, 1980

#  - A particular model formulation that was originally designed for evaluating diagnostic tests in the absence of a gold standard

#- Not originally/necessarily Bayesian - implemented using Maximum Likelihood 

# - Evaluating an imperfect test against another imperfect test; is a bit like pulling a rabbit out of a hat

# If we don't know the true disease status, how can we estimate sensitivity or specificity for either test?

# https://www.youtube.com/watch?v=z6devQmW2xE&ab_channel=PolychronisKostoulas#



# We will use the data/observations from the manuscript published back in 1980.

## Hui-Walter (1980) dataset Table 1
pop_1 = matrix(nrow=3,ncol=3)
rownames(pop_1) = c("Mantoux_Test_Pos", "Mantoux_Test_Neg", "Total")
colnames(pop_1) = c("Tine_Test_Pos", "Tine_Test_Neg", "Total")

pop_1[1,1] = 14
pop_1[1,2] = 4
pop_1[2,1] = 9
pop_1[2,2] = 528
#Total rows and columns
pop_1[1,3] = pop_1[1,1] + pop_1[1,2]
pop_1[2,3] = pop_1[2,1] + pop_1[2,2]
pop_1[3,1] = pop_1[1,1] + pop_1[2,1]
pop_1[3,2] = pop_1[1,2] + pop_1[2,2]
N_1 = sum(pop_1[1,1] + pop_1[1,2] + pop_1[2,1] + pop_1[2,2])
pop_1[3,3] = N_1
pop_1

## Now let's do pop_2
pop_2 = matrix(nrow=3,ncol=3)
rownames(pop_2) = c("Mantoux_Test_Pos", "Mantoux_Test_Neg", "Total")
colnames(pop_2) = c("Tine_Test_Pos", "Tine_Test_Neg", "Total")

pop_2[1,1] = 887
pop_2[1,2] = 31
pop_2[2,1] = 37
pop_2[2,2] = 367
#Total rows and columns
pop_2[1,3] = pop_2[1,1] + pop_2[1,2]
pop_2[2,3] = pop_2[2,1] + pop_2[2,2]
pop_2[3,1] = pop_2[1,1] + pop_2[2,1]
pop_2[3,2] = pop_2[1,2] + pop_2[2,2]
N_2 = sum(pop_2[1,1] + pop_2[1,2] + pop_2[2,1] + pop_2[2,2])
pop_2[3,3] = N_2
pop_2

## Exercise

# Assuming Mantoux test as a gold standard, estimate and save the sensitivity and specificity of tine test in both populations?


## Solution

# 1st population
(sensitivity_1 <- pop_1[1,1] / (pop_1[1,3]))
(specificity_1 <- pop_1[2,2] / (pop_1[2,3]))

#2nd population
(sensitivity_2 <- pop_2[1,1] / (pop_2[1,3]))
(specificity_2 <- pop_2[2,2] / (pop_2[2,3]))

## Hui-Walter model

#  - A particular model formulation that was originally designed for evaluating diagnostic tests in the absence of a gold standard

# - Also known as the two_test - two_population setting/paradigm



## Model Specification ('hw_definition')

hw_definition <- c("model{
  Population_1 ~ dmulti(prob_1, N_1)
  Population_2 ~ dmulti(prob_2, N_2)
  
  #Population_1
  
  # Test1+ Test2+
	prob_1[1] <- (prev[1] * ((se[1])*(se[2]))) + ((1-prev[1]) * ((1-sp[1])*(1-sp[2])))
  
  # Test1+ Test2-
	prob_1[2] <- (prev[1] * ((se[1])*(1-se[2]))) + ((1-prev[1]) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob_1[3] <- (prev[1] * ((1-se[1])*(se[2]))) + ((1-prev[1]) * ((sp[1])*(1-sp[2])))

  # Test1- Test2-
	prob_1[4] <- (prev[1] * ((1-se[1])*(1-se[2]))) + ((1-prev[1]) * ((sp[1])*(sp[2])))
	
	#Population_2
  
  # Test1+ Test2+
	prob_2[1] <- (prev[2] * ((se[1])*(se[2]))) + ((1-prev[2]) * ((1-sp[1])*(1-sp[2])))
  
  # Test1+ Test2-
	prob_2[2] <- (prev[2] * ((se[1])*(1-se[2]))) + ((1-prev[2]) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob_2[3] <- (prev[2] * ((1-se[1])*(se[2]))) + ((1-prev[2]) * ((sp[1])*(1-sp[2])))

  # Test1- Test2-
	prob_2[4] <- (prev[2] * ((1-se[1])*(1-se[2]))) + ((1-prev[2]) * ((sp[1])*(sp[2])))

  prev[1] ~ dbeta(1, 1)
  prev[2] ~ dbeta(1, 1)
  
  se[1] ~ dbeta(1, 1)I(1-sp[1], )
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)I(1-sp[2], )
  sp[2] ~ dbeta(1, 1)

  #data# Population_1, Population_2, N_1, N_2
  #monitor# prev, prob_1, prob_2, se, sp
  #inits# prev, se, sp
  }
  ")
library('runjags')

Population_1 <- as.numeric(pop_1[1:2,1:2])
Population_2 <- as.numeric(pop_2[1:2,1:2])


prev <- list(chain1=c(0.05,0.99), chain2=c(0.95,0.05))
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))

results <- run.jags(hw_definition, n.chains=2)


# Remember to check convergence and effective sample size!

plot(results)

pt <- plot(results)
pt$`prev[1].plot1`
pt$`prev[1].plot3`

print(pt[["prev[1].plot1"]])

print(pt[["se[1].plot1"]])
print(pt[["sp[1].plot1"]])
print(pt[["sp[1].plot3"]])

summary(results)

## Exercise 1 

#  Run the `hw_definition` model under the following different scenarios 
# and interpret the results in each case.

# 1. Change the priors for *Se[1]* and *Sp[1]* and try Beta(5,1).
curve(dbeta(x,5,1))  # Plots the distribution of Beta(2,1)


# Solution
## Model Specification ('hw_definition')

hw_definition_1 <- c("model{
  Population_1 ~ dmulti(prob_1, N_1)
  Population_2 ~ dmulti(prob_2, N_2)
  
  #Population_1
  
  # Test1+ Test2+
	prob_1[1] <- (prev[1] * ((se[1])*(se[2]))) + ((1-prev[1]) * ((1-sp[1])*(1-sp[2])))
  
  # Test1+ Test2-
	prob_1[2] <- (prev[1] * ((se[1])*(1-se[2]))) + ((1-prev[1]) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob_1[3] <- (prev[1] * ((1-se[1])*(se[2]))) + ((1-prev[1]) * ((sp[1])*(1-sp[2])))

  # Test1- Test2-
	prob_1[4] <- (prev[1] * ((1-se[1])*(1-se[2]))) + ((1-prev[1]) * ((sp[1])*(sp[2])))
	
	#Population_2
  
  # Test1+ Test2+
	prob_2[1] <- (prev[2] * ((se[1])*(se[2]))) + ((1-prev[2]) * ((1-sp[1])*(1-sp[2])))
  
  # Test1+ Test2-
	prob_2[2] <- (prev[2] * ((se[1])*(1-se[2]))) + ((1-prev[2]) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob_2[3] <- (prev[2] * ((1-se[1])*(se[2]))) + ((1-prev[2]) * ((sp[1])*(1-sp[2])))

  # Test1- Test2-
	prob_2[4] <- (prev[2] * ((1-se[1])*(1-se[2]))) + ((1-prev[2]) * ((sp[1])*(sp[2])))

  prev[1] ~ dbeta(1, 1)
  prev[2] ~ dbeta(1, 1)
  
  se[1] ~ dbeta(5, 1)I(1-sp[1], )
  sp[1] ~ dbeta(5, 1)
  se[2] ~ dbeta(1, 1)I(1-sp[2], )
  sp[2] ~ dbeta(1, 1)

  #data# Population_1, Population_2, N_1, N_2
  #monitor# prev, prob_1, prob_2, se, sp
  #inits# prev, se, sp
  }
  ")

Population_1 <- as.numeric(pop_1[1:2,1:2])
Population_2 <- as.numeric(pop_2[1:2,1:2])


prev <- list(chain1=c(0.05,0.99), chain2=c(0.95,0.05))
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))

results_1 <- run.jags(hw_definition_1, n.chains=2)
# Remember to check convergence and effective sample size!
plot(results_2)  # Click backwards to view all plots

pt_1 <- plot(results_1)
pt_1$`prev[1].plot1`
pt_1$`prev[1].plot3`

print(pt_1[["prev[1].plot1"]])

print(pt_1[["se[1].plot1"]])
print(pt_1[["sp[1].plot1"]])
print(pt_1[["sp[1].plot3"]])

summary(results_1)

# 2. Remove the `I(1-sp[1], )` and 'I(1-sp[2])' from the model and run it again. What happens now?

# Solution
## Model Specification ('hw_definition')

hw_definition_2 <- c("model{
  Population_1 ~ dmulti(prob_1, N_1)
  Population_2 ~ dmulti(prob_2, N_2)
  
  #Population_1
  
  # Test1+ Test2+
	prob_1[1] <- (prev[1] * ((se[1])*(se[2]))) + ((1-prev[1]) * ((1-sp[1])*(1-sp[2])))
  
  # Test1+ Test2-
	prob_1[2] <- (prev[1] * ((se[1])*(1-se[2]))) + ((1-prev[1]) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob_1[3] <- (prev[1] * ((1-se[1])*(se[2]))) + ((1-prev[1]) * ((sp[1])*(1-sp[2])))

  # Test1- Test2-
	prob_1[4] <- (prev[1] * ((1-se[1])*(1-se[2]))) + ((1-prev[1]) * ((sp[1])*(sp[2])))
	
	#Population_2
  
  # Test1+ Test2+
	prob_2[1] <- (prev[2] * ((se[1])*(se[2]))) + ((1-prev[2]) * ((1-sp[1])*(1-sp[2])))
  
  # Test1+ Test2-
	prob_2[2] <- (prev[2] * ((se[1])*(1-se[2]))) + ((1-prev[2]) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob_2[3] <- (prev[2] * ((1-se[1])*(se[2]))) + ((1-prev[2]) * ((sp[1])*(1-sp[2])))

  # Test1- Test2-
	prob_2[4] <- (prev[2] * ((1-se[1])*(1-se[2]))) + ((1-prev[2]) * ((sp[1])*(sp[2])))

  prev[1] ~ dbeta(1, 1)
  prev[2] ~ dbeta(1, 1)
  
  se[1] ~ dbeta(1, 1)
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)
  sp[2] ~ dbeta(1, 1)

  #data# Population_1, Population_2, N_1, N_2
  #monitor# prev, prob_1, prob_2, se, sp
  #inits# prev, se, sp
  }
  ")

Population_1 <- as.numeric(pop_1[1:2,1:2])
Population_2 <- as.numeric(pop_2[1:2,1:2])


prev <- list(chain1=c(0.05,0.99), chain2=c(0.95,0.05))
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))

results_2 <- run.jags(hw_definition_2, n.chains=2)
# Remember to check convergence and effective sample size!
plot(results_2)  # Click backwards to view all plots

pt_2 <- plot(results_2)
pt_2$`prev[1].plot1`
pt_2$`prev[1].plot3`

print(pt_2[["prev[1].plot1"]])

print(pt_2[["se[1].plot1"]])
print(pt_2[["sp[1].plot1"]])
print(pt_2[["sp[1].plot3"]])

print(pt_2[["se[2].plot1"]])

summary(results_2)

# 3. Try to run the model with different initial values, that explore the whole parameter space and removing `T(1-sp[1],)`.
# For example try it with:
# se <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))
# sp <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))

# Solution
## Model Specification ('hw_definition')
# The model is the same like as in the above scenario (hw_definition_2)

Population_1 <- as.numeric(pop_1[1:2,1:2])
Population_2 <- as.numeric(pop_2[1:2,1:2])


prev <- list(chain1=c(0.05,0.99), chain2=c(0.95,0.05))
se <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))
sp <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))

results_3 <- run.jags(hw_definition_2, n.chains=2, sample = 100000)
# Remember to check convergence and effective sample size!
plot(results_3)  # Click backwards to view all plots

pt_3 <- plot(results_3)
pt_3$`prev[1].plot1`
pt_3$`prev[1].plot3`

print(pt_3[["prev[1].plot1"]])

print(pt_3[["se[1].plot1"]])
print(pt_3[["sp[1].plot1"]])
print(pt_3[["sp[1].plot3"]])
print(pt_3[["se[2].plot1"]])

summary(results_3)

# The model fails to converge to one solution, but reports two solution that are complementary

# 4. Run the model with only 1 population (either pop_1 or pop_2). What happens then?
# We'll remove population 2


# Solution
## Model Specification ('hw_definition')

hw_definition_4 <- c("model{
  Population_1 ~ dmulti(prob_1, N_1)

  #Population_1
  
  # Test1+ Test2+
	prob_1[1] <- (prev[1] * ((se[1])*(se[2]))) + ((1-prev[1]) * ((1-sp[1])*(1-sp[2])))
  
  # Test1+ Test2-
	prob_1[2] <- (prev[1] * ((se[1])*(1-se[2]))) + ((1-prev[1]) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob_1[3] <- (prev[1] * ((1-se[1])*(se[2]))) + ((1-prev[1]) * ((sp[1])*(1-sp[2])))

  # Test1- Test2-
	prob_1[4] <- (prev[1] * ((1-se[1])*(1-se[2]))) + ((1-prev[1]) * ((sp[1])*(sp[2])))
	

  prev[1] ~ dbeta(1, 1)
  se[1] ~ dbeta(1, 1)I(1-sp[1],)
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)I(1-sp[2],)
  sp[2] ~ dbeta(1, 1)

  #data# Population_1, N_1
  #monitor# prev, prob_1,  se, sp
  #inits# prev, se, sp
  }
  ")

Population_1 <- as.numeric(pop_1[1:2,1:2])


prev <- list(chain1=c(0.05), chain2=c(0.95))
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))

results_4 <- run.jags(hw_definition_4, n.chains=2)
# Remember to check convergence and effective sample size!
plot(results_4)  # Click backwards to view all plots

pt_4 <- plot(results_4)
pt_4$`prev.plot1`
pt_4$`prev.plot3`

print(pt_4[["prev.plot1"]])

print(pt_4[["se[1].plot1"]])
print(pt_4[["sp[1].plot1"]])
print(pt_4[["sp[1].plot3"]])
print(pt_4[["se[2].plot3"]])
print(pt_4[["se[2].plot4"]])
print(pt_4[["se[2].plot5"]])

summary(results_4)
# How do the results look? # Compare the results with the scenarios above? 
# Try also by removing population 2




#################################

# References for Reading

# Estimating the Error Rates of Diagnostic Tests, S. L. Hui and S. D. Walter, 1980 (https://doi.org/10.2307/2530508)
# Bayesian estimation of disease prevalence and the parameters of diagnostic tests in the absence of a gold standard, Joseph et al, 1995 (https://doi.org/10.1093/oxfordjournals.aje.a117428)
# STARD-BLCM: Standards for the Reporting of Diagnostic accuracy studies that use Bayesian Latent Class Models, Kostoulas et al, 2017 (https://doi.org/10.1016/j.prevetmed.2017.01.006)
# Diagnostic Accuracy Estimates for COVID-19 Real-Time Polymerase Chain Reaction and Lateral Flow Immunoassay Tests With Bayesian Latent-Class Models, Kostoulas et al, 2021 (https://doi.org/10.1093/aje/kwab093)
# |tPRiors |: a tool for prior elicitation and obtaining posterior distributions of true disease prevalence, Pateras - Kostoulas 2022 (https://doi.org/10.1186/s12874-022-01557-1)










