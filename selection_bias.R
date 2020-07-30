#-----------------------------------------------#
# Data example for Quality over Quantity article
#-----------------------------------------------#

#---------------------#
# Dependencies
#---------------------#
rm(list = ls())
library(tidyverse)
library(bindata)

# Set working directory below:
# setwd("")

#-----------------------------#
# Establish joint distribution
#-----------------------------#

# Parameters of joint distribution
p_symptoms = 0.5
p_infected = 0.35

# Check range of correlations via odds ratios
bincorr <- function(OR, p1, p2) {    #from odds ratio to binary correlation
  if (OR==1) p11=p2-p2+p1*p2 else {
    p11_1=p2-(1/2/(1-OR)*(1-p1+OR*p1+p2-OR*p2-
                            sqrt((-1+p1-OR*p1-p2+OR*p2)^2-4*(1-OR)*(p2-p1*p2))))
    p11_2=p2-(1/2/(1-OR)*(1-p1+OR*p1+p2-OR*p2-
                            sqrt((-1+p1-OR*p1-p2+OR*p2)^2)-4*(1-OR)*(p2-p1*p2)))
    if (p11_1>0 && p11_1<=p1 && p11_1<p2) p11=p11_1 else p11=p11_2
  }
  bincorr=(p11-p1*p2)/sqrt(p1*(1-p1)*p2*(1-p2))
  return(bincorr)
}

sapply(c(0,0.5,1,10,100,1000),function(x) bincorr(x,p_symptoms,p_infected))

# Use probabilities above, and a correlation between symptomatic and infection to generate population of patients
size = 1e6
rho = 0.7
population = rmvbin(size, c(p_symptoms,p_infected), bincorr = (1-rho) * diag(2) + rho)

# Examine truly random data model for consistency
skimr::skim( population[,1] )
skimr::skim( population[,2] )

#---------------------------#
# Trouble with selection bias
#---------------------------#

# Examine non-random selection model for consistency
subset_symptomatic = population[ population[,1] == 1 , ]

# Different sample sizes from non-randomly sampled data
samp_size = c(10 , 100 , 1000 , 10000 , 100000)

for (s in 1:length(samp_size)) {
  print( skimr::skim( subset_symptomatic[ sample( nrow(subset_symptomatic) , samp_size[s] ) , 2 ] ) )
}

#----------------------------#
# Compare with random sampling
#----------------------------#

# Different sample sizes from randomly sampled data
for (s in 1:length(samp_size)) {
  print( skimr::skim( population[ sample( nrow(population) , samp_size[s] ) , 2 ] ) )
}

#-------------------------------------#
# Plot biased and randomly sampled data
#-------------------------------------#

# Write a function to graphically display effects of selective sampling bias
prop_function = function(rho) {
  
  # Use probabilities above and correlation between symptomatic and infection to generate population of patients
  size = 1e6
  p_symptoms = 0.5
  p_infected = 0.35
  population = rmvbin(size, c(p_symptoms,p_infected), bincorr = (1-rho) * diag(2) + rho)

  # Sample size
  nsize = 5e3

  # Sample from symptomatic population
  subset_symptomatic = population[ population[,1] == 1 , ]
  select_samp = subset_symptomatic[ sample( nrow(subset_symptomatic) , nsize ) , 2 ]

  # Sample randomly
  random_samp = population[ sample( nrow(population) , nsize ) , 2 ]

  # Accumulating data
  index = c(1:nsize)
  random = cumsum(random_samp) / index
  symptomatic = cumsum(select_samp) / index

  # Plot objects above
  plot( index , random , type = "l" , ylim = c(0,1) , 
        xlab = "Tested" , ylab = "Infected Proportion" , cex.axis = 2 , cex.lab = 2.5 )
  lines( symptomatic , col = "light blue" )
  abline(h = p_infected , lty = 2)
  title( main = paste0("Correlation: " , round(rho,2)) , cex.main = 4 )
}

# Graphic for what happens as sampling grows under different correlations between symptomatic and infected patients
png("rand_sympt_sample_bias.png" , width = 1500 , height = 1500)
par( mfrow = c(3,3) , mar = c(6,6,4,4) )
sapply( seq(0 , 0.7 , by = 0.1) , prop_function )
dev.off()

pdf("rand_sympt_sample_bias.pdf")
par( mfrow = c(3,3) )
sapply( seq(0 , 0.7 , by = 0.1) , prop_function )
dev.off()
