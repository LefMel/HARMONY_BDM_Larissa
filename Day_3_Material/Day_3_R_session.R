
library("runjags")

# Sensitivity - Specificity estimation with and without a gold standard

## Hui-Walter paradigm/model (1980)
  
#  - A particular model formulation that was originally designed for evaluating diagnostic tests in the absence of a gold standard

#- Not originally/necessarily Bayesian - implemented using Maximum Likelihood 

# - Evaluating an imperfect test against another imperfect test; is a bit like pulling a rabbit out of a hat

# If we don't know the true disease status, how can we estimate sensitivity or specificity for either test?

#https://www.youtube.com/watch?v=z6devQmW2xE&ab_channel=PolychronisKostoulas#

## Hui-Walter paradigm (1980)

# * Hui-Walter models implementation to be further discussed in the next session.

# * But we will use the data/observations from the manuscript published back in 1980.

## Hui-Walter (1980) dataset

## Encode the Table_1 data in RStudio


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

![](figs/hui.walter.pdf)


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
  
  se[1] ~ dbeta(1, 1)T(1-sp[1], )
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)T(1-sp[2], )
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
print(pt[["prev[1].plot1"]])

print(pt[["se[1].plot1"]])
print(pt[["sp[1].plot1"]])
print(pt[["sp[1].plot3"]])

summary(results)

  ## Exercise 1 
  
#  Run the `hw_definition` model under the following different scenarios and interpret the results in each case.

# 1. Change the priors for *Se[1]* and *Sp[1]* and try Beta(2,1).

# 2. Remove the `T(1-sp[1], )` from the model and run it again. What happens now?
  
# 3. Try to run the model with different initial values and removing `T(1-sp[1],)`.
#For example try it with:
se <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))
sp <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))

# 4. Run the model with only 1 population (either pop_1 or pop_2). What happens then?
  
  
## Event Closure
  
#  Group photo







