#MLPE models for landscape genetic analysis of Circuitscape resistance in Emel et al. (2019)
#Sarah Emel

library(MuMIn)
library(lme4)

##this is an example. run separate models for each study area, genetic distance, and set of landscape variables##

#genetic distance data
dps <- read.csv(file = "rvar1dps.csv", na.strings = "na")
dps <- dps[lower.tri(dps)]

#populations
pops<-read.csv(file="rv1pops.csv",header=T)
pop1 <- factor(pops$pop1)
pop2 <- factor(pops$pop2)

#circuitscape cost distance
cs <- read.csv(file = "CSout/rv1cs.csv")
cs <- cs[,-1] #remove column of row names
cs_scaled <- data.frame(lapply(cs, FUN = scale))

##select response variable##
gdist <- dps

#create dataframe
data <- cbind(pop1, pop2, gdist, cs_scaled)

##formula part 1
Zl <- lapply(c("pop1","pop2"), function(nm) Matrix:::fac2sparse(data[[nm]],"d", drop=FALSE))
ZZ <- Reduce("+", Zl[-1], Zl[[1]])

##formula part 2

#MLPE functions
#use MLPE for model fitting, examining coefficients, CIs
#use MLPEnoREML for model comparison
MLPE <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = TRUE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

MLPEnoREML <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = FALSE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

#do model validation with MLPE function
var <- data$rv1canor #specify variable
fit <- MLPE(gdist ~ var + (1|pop1), data=data) 
hist(residuals(fit)) #look for roughly normal distribution
plot(fit) #look for cone shape

#heading for the output file
header <- noquote(paste("AICc", "BIC", "R2m", "R2c", "Var",sep=" ")) 
write.table(header, file = "rv1csdps.txt", append=T, row.names=F, col.names=F, quote = F) 

#do model comparison with MLPEnoREML function
for (i in (4:ncol(data))) {
  fit <- MLPEnoREML(gdist ~ data[, i] + (1|pop1), data=data)
  anova(fit)
  b <- BIC(fit)
  a <- AICc(fit, k = length(unique(pop1))) #this make k = # pops, not # pairs
  r <- r.squaredGLMM(fit)
  rm <- r[1]
  rc <- r[2]
  out<- noquote(paste(a, b, rm, rc, colnames(data)[i], sep =" "))
  write.table(out ,file = "rv1csdps.txt", append=T,row.names=F,col.names=F, quote = F) 
}

#get confidence intervals for predictors
header2 <- noquote(paste("Var","Est","2.5%", "97.5%", sep=" ")) 
write.table(header2, file = "rv1csdps_ci.txt", append=T, row.names=F, col.names=F, quote = F)

for (i in (4:ncol(data))) { #may need to adjust code for multivariate models
  fit <- MLPE(gdist ~ data[, i] + (1|pop1), data=data)
  est <- coef(summary(fit))
  ci <- confint(fit, level = 0.95, method = "Wald")
  out<- noquote(paste(colnames(data)[i], est[2,1], ci[4,1], ci[4,2], sep =" "))
  write.table(out ,file = "rv1csdps_ci.txt", append=T,row.names=F,col.names=F, quote = F) 
}