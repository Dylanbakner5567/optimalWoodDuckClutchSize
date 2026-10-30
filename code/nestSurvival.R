
library(jagsUI)

#setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#set your working directory to the project root before running this script
#the folder should contain: data/ (with jags.data2.rds) and code/ (for model output)
#setwd("path/to/project/root")

#load jags data
jags.data2 <- readRDS("data/jags.data2.rds")

#write jags model file to the working directory
cat(file="model2.txt", "
model {

#likelihood~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#nest survival model

    for(i in 1:nClutches){
      for(t in (first[i]+1):last[i]){
        y[i,t] ~ dbern(p[i,t])
        logit(q[i,t]) <- mean.b0 + beta1*cs[i,t] + beta2*cs2[i,t]
        p[i,t] <- (y[i,t-1])*q[i,t]
      }
    }

#priors~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  b0 ~ dunif(0,1)
  mean.b0 <- logit(b0)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)

}
")

inits <- function(){list()}

#parameters monitored
parameters <- c("b0", "mean.b0", "beta1", "beta2")

#mcmc settings
ni <-80000; nb <- 30000; nc <- 4; nt <- 10

#call jags from r
m1 <- jags(jags.data2, inits, parameters, "model2.txt", n.chains = nc,
          n.thin = nt, n.iter = ni, n.burnin = nb, parallel = F)
summary(m1)

