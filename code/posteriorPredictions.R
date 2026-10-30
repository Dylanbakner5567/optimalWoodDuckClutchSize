
library(jagsUI)
library(dplyr)

#setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#set your working directory to the project root before running this script
#the folder should contain: data/ (with jags.data2.rds and resultData.RData)
#setwd("path/to/project/root")

#load data and model output
jags.data2 <- readRDS("data/jags.data2.rds")
load("data/resultData.RData")
m  <- resultData$m
m1 <- resultData$m1

#daily growth rates~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#build day and lambda matrices aligned to each nest's active period
day <- matrix(NA, ncol = ncol(m$mean$lambda), nrow = nrow(m$mean$lambda))
dailyGrowth <- matrix(NA, ncol = ncol(m$mean$lambda), nrow = nrow(m$mean$lambda))
for (i in 1:nrow(day)) {
  start <- resultData$nestInit[i]
  end   <- resultData$julChecked[i]
  len   <- end - start + 1
  day[i, start:end]         <- 1:len
  dailyGrowth[i, start:end] <- m$mean$lambda[i, start:end]
}

#reshape to long format
dayV    <- as.vector(unlist(day))
growthV <- as.vector(unlist(dailyGrowth))
dailyGrowthPlot <- data.frame(day = dayV, growth = growthV)
dailyGrowthPlot <- na.omit(dailyGrowthPlot)

#quantiles across all nests and all days
dailyGrowthPlot %>% summarise(
  qLow  = quantile(growth, 0.025, na.rm = TRUE),
  qMed  = quantile(growth, 0.5,   na.rm = TRUE),
  qHigh = quantile(growth, 0.975, na.rm = TRUE)
)

#quantiles during the laying period (days 1-14)
layGrowthVals <- dailyGrowthPlot %>%
  filter(day %in% 1:14) %>%
  summarise(
    qLow  = quantile(growth, 0.025, na.rm = TRUE),
    qMed  = quantile(growth, 0.5,   na.rm = TRUE),
    qHigh = quantile(growth, 0.975, na.rm = TRUE)
  )
layGrowthVals

#quantiles during the incubation period (days 15-45)
incGrowthValues <- dailyGrowthPlot %>%
  filter(day %in% 15:45) %>%
  summarise(
    qLow  = quantile(growth, 0.025, na.rm = TRUE),
    qMed  = quantile(growth, 0.5,   na.rm = TRUE),
    qHigh = quantile(growth, 0.975, na.rm = TRUE)
  )
incGrowthValues

#daily growth rates for successful nests~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#rebuild matrices (same structure as above)
day2 <- matrix(NA, ncol = ncol(m$mean$lambda), nrow = nrow(m$mean$lambda))
dailyGrowth2 <- matrix(NA, ncol = ncol(m$mean$lambda), nrow = nrow(m$mean$lambda))
for (i in 1:nrow(day2)) {
  start <- resultData$nestInit[i]
  end   <- resultData$julChecked[i]
  len   <- end - start + 1
  day2[i, start:end]         <- 1:len
  dailyGrowth2[i, start:end] <- m$mean$lambda[i, start:end]
}

#retain only nests where all non-NA survival values equal 1
keep <- which(apply(jags.data2$y, 1, function(x) any(!is.na(x)) && all(x[!is.na(x)] == 1)))
dailyGrowth2 <- dailyGrowth2[keep, ]

#reshape to long format
dayV2    <- as.vector(unlist(day2[keep, ]))
growthV2 <- as.vector(unlist(dailyGrowth2))
dailyGrowthPlot2 <- data.frame(day = dayV2, growth = growthV2)
dailyGrowthPlot2 <- na.omit(dailyGrowthPlot2)

#quantiles across all days for successful nests
dailyGrowthPlot2 %>% summarise(
  qLow  = quantile(growth, 0.025, na.rm = TRUE),
  qMed  = quantile(growth, 0.5,   na.rm = TRUE),
  qHigh = quantile(growth, 0.975, na.rm = TRUE)
)

#quantiles during laying for successful nests
layGrowthVals <- dailyGrowthPlot2 %>%
  filter(day %in% 1:14) %>%
  summarise(
    qLow  = quantile(growth, 0.025, na.rm = TRUE),
    qMed  = quantile(growth, 0.5,   na.rm = TRUE),
    qHigh = quantile(growth, 0.975, na.rm = TRUE)
  )
layGrowthVals

#quantiles during incubation for successful nests
incGrowthValues <- dailyGrowthPlot2 %>%
  filter(day %in% 15:45) %>%
  summarise(
    qLow  = quantile(growth, 0.025, na.rm = TRUE),
    qMed  = quantile(growth, 0.5,   na.rm = TRUE),
    qHigh = quantile(growth, 0.975, na.rm = TRUE)
  )
incGrowthValues

#nest survival - clutch size effect~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#posterior means and credible intervals for clutch size coefficients
m1$mean$beta1; m1$q2.5$beta1; m1$q97.5$beta1
m1$mean$beta2; m1$q2.5$beta2; m1$q97.5$beta2

#scaling constants derived from observed clutch sizes
muCS <- mean(resultData$clutchMat, na.rm = TRUE)
sdCS <- sd(resultData$clutchMat,   na.rm = TRUE)

#standardize clutch sizes 1-30 for prediction
nScale <- vector()
for (i in 1:30) {
  nScale[i] <- (i - muCS) / sdCS
}
nScale2 <- nScale^2

#predict daily survival probability across clutch sizes using posterior draws
ns <- matrix(data = NA, ncol = length(m1$sims.list$b0), nrow = 30)
for (d in 1:nrow(ns)) {
  for (i in 1:length(ns[1, ])) {
    ns[d, i] <- m1$sims.list$mean.b0[i] +
      m1$sims.list$beta1[i] * nScale[d] +
      m1$sims.list$beta2[i] * nScale2[d]
  }
}
ns <- plogis(ns)

#summarise predictions across posterior draws
quants <- t(apply(ns, 1, quantile, probs = c(0.025, 0.5, 0.975)))
predNS <- data.frame(cs     = 1:nrow(ns),
                     lower  = quants[, 1],
                     median = quants[, 2],
                     upper  = quants[, 3])

#clutch size that maximizes daily survival probability
predNS[which.max(predNS$median), ]

#scenario predictions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#q25 (loss) and q75 (gain) growth rates by day for scenario inputs
lossValsL <- data.frame(dailyGrowthPlot %>%
  filter(day %in% 1:14) %>%
  group_by(day) %>%
  summarise(qMed = quantile(growth, 0.25, na.rm = TRUE)))

lossValsI <- data.frame(dailyGrowthPlot %>%
  filter(day %in% 15:45) %>%
  group_by(day) %>%
  summarise(qMed = quantile(growth, 0.25, na.rm = TRUE)))

gainValsL <- data.frame(dailyGrowthPlot %>%
  filter(day %in% 1:14) %>%
  group_by(day) %>%
  summarise(qMed = quantile(growth, 0.75, na.rm = TRUE)))

gainValsI <- data.frame(dailyGrowthPlot %>%
  filter(day %in% 15:45) %>%
  group_by(day) %>%
  summarise(qMed = quantile(growth, 0.75, na.rm = TRUE)))

#function to predict daily survival probability for a given clutch trajectory
makePred <- function(layVals, incVals) {

  #laying stage (days 1-14): accumulate clutch size day by day
  i <- 0:13
  n <- i + (1 + cumsum(layVals[1:14]))

  #standardize
  nScale  <- (n - muCS) / sdCS
  nScale2 <- nScale^2

  #predict using posterior draws
  ns <- sapply(1:length(m1$sims.list$mean.b0), function(s) {
    m1$sims.list$mean.b0[s] +
      m1$sims.list$beta1[s] * nScale +
      m1$sims.list$beta2[s] * nScale2
  })
  ns <- plogis(ns)
  q1    <- t(apply(ns, 1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE))
  pred1 <- data.frame(day = 1:14, lower = q1[, 1], median = q1[, 2], upper = q1[, 3], clutch = n)

  #incubation stage (days 15-45): continue accumulating from end of laying
  n1       <- n[14] + cumsum(incVals[1:31])
  nScale1  <- (n1 - muCS) / sdCS
  nScale21 <- nScale1^2

  ns1 <- sapply(1:length(m1$sims.list$mean.b0), function(s) {
    m1$sims.list$mean.b0[s] +
      m1$sims.list$beta1[s] * nScale1 +
      m1$sims.list$beta2[s] * nScale21
  })
  ns1 <- plogis(ns1)
  q2    <- t(apply(ns1, 1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE))
  pred2 <- data.frame(day = 15:45, lower = q2[, 1], median = q2[, 2], upper = q2[, 3], clutch = n1)

  rbind(pred1, pred2)
}

#generate predictions for four scenarios
predConstantDep     <- makePred(lossValsL[, 2], lossValsI[, 2])
predConstantPara    <- makePred(gainValsL[, 2], gainValsI[, 2])
predHighDepHighPara <- makePred((lossValsL[, 2] + gainValsL[, 2]), (lossValsI[, 2] + gainValsI[, 2]))
predNoDepNoPara     <- makePred(rep(0, 14), rep(0, 31))

#scenario nest survival probabilities~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#cumulative survival probability during laying (days 1-14)
prod(predNoDepNoPara$median[1:14]);     prod(predNoDepNoPara$lower[1:14]);     prod(predNoDepNoPara$upper[1:14])
prod(predConstantPara$median[1:14]);    prod(predConstantPara$lower[1:14]);    prod(predConstantPara$upper[1:14])
prod(predConstantDep$median[1:14]);     prod(predConstantDep$lower[1:14]);     prod(predConstantDep$upper[1:14])
prod(predHighDepHighPara$median[1:14]); prod(predHighDepHighPara$lower[1:14]); prod(predHighDepHighPara$upper[1:14])

#cumulative survival probability during incubation (days 15-45)
prod(predNoDepNoPara$median[15:45]);     prod(predNoDepNoPara$lower[15:45]);     prod(predNoDepNoPara$upper[15:45])
prod(predConstantPara$median[15:45]);    prod(predConstantPara$lower[15:45]);    prod(predConstantPara$upper[15:45])
prod(predConstantDep$median[15:45]);     prod(predConstantDep$lower[15:45]);     prod(predConstantDep$upper[15:45])
prod(predHighDepHighPara$median[15:45]); prod(predHighDepHighPara$lower[15:45]); prod(predHighDepHighPara$upper[15:45])

#overall cumulative survival probability (days 1-45)
prod(predNoDepNoPara$median);     prod(predNoDepNoPara$lower);     prod(predNoDepNoPara$upper)
prod(predConstantPara$median);    prod(predConstantPara$lower);    prod(predConstantPara$upper)
prod(predConstantDep$median);     prod(predConstantDep$lower);     prod(predConstantDep$upper)
prod(predHighDepHighPara$median); prod(predHighDepHighPara$lower); prod(predHighDepHighPara$upper)
