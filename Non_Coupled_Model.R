#code for the Non-coupled Model from Section 4 of the main text
#set your working directory to source file location
rm(list=ls())
library(questionr)
library(dplyr)
library(nimble)
library(coda)
library(tidyr)
library(padr)
library(stringr)
library(lsr)

memory.limit(size=1E10)

#load in the data
load("Quebec_hospital_study_data.RData")
mT <- length(admissions[1,])
N <- length(admissions[,1])
lpsi <- log(admissions+1)
mlpsi <- mean(lpsi)

#Nimble constants and data
hospitalConsts <- list(mT=length(admissions[1,]),N=length(admissions[,1]),lpsi=lpsi,mlpsi=mlpsi,
                       standard_cap=standard_cap,
                       mobility_matrix=mobility_matrix,mmobility_matrix=mean(mobility_matrix),
                       new_variant=new_variant)

hospitalData <- list(y=admissions,constraint_data1=1,
                     constraint_data2=1,constraint_data3=1,
                     constraint_data4=1)


#Nimble model code, note that not all variables are named the same as in the main text
hospitalCode <- nimbleCode({
  
  #likelihood
  indi[1:7] <- c(0,1,1,1,1,1,1)
  for(i in 1:N){
    for(t in 2:mT){
      
      #data likelihood
      y[i,t] ~ dnegbin(prob=r[i,t,S[i,t]]/(r[i,t,S[i,t]]+lambda[i,t,S[i,t]]*indi[S[i,t]])-1e-10*(1-indi[S[i,t]]),
                       size=r[i,t,S[i,t]])
      
      lambda[i,t,1:7] <- c(1,lambda_en[i,t],lambda_en[i,t],lambda_ep[i,t],lambda_ep[i,t],lambda_ep[i,t],lambda_ep[i,t])
      r[i,t,1:7] <- c(1,r0,r0,r1,r1,r1,r1)
      lambda_ep[i,t] <- exp(b1[i]+
                              rho1*(lpsi[i,t-1]-mlpsi))
      lambda_en[i,t] <-  exp(b0[i]+
                               rho0*(lpsi[i,t-1]-mlpsi))
    }
  }
  
  #Markov chain
  uniff[1:7] <- c(1/3,1/6,1/6,1/12,1/12,1/12,1/12)
  epi_indi[1:7] <- c(0,0,0,1,1,1,1)
  for(i in 1:N){
    S[i,1] ~ dcat(uniff[1:7])
    for(t in 2:mT){
      
      #calc transition matrix
      tm[i,t,1,1:7] <- c(1-p12[i,t],p12[i,t],0,0,0,0,0)
      tm[i,t,2,1:7] <- c(0,0,1,0,0,0,0)
      tm[i,t,3,1:7] <- c(p21[i,t],0,p22[i,t],p23[i,t],0,0,0)
      tm[i,t,4,1:7] <- c(0,0,0,0,1,0,0)
      tm[i,t,5,1:7] <- c(0,0,0,0,0,1,0)
      tm[i,t,6,1:7] <- c(0,0,0,0,0,0,1)
      tm[i,t,7,1:7] <- c(0,1-p33[i,t],0,0,0,0,p33[i,t])
      logit(p12[i,t]) <- alpha[1]+alpha[5]*(standard_cap[i]-mean(standard_cap[1:N]))+
        alpha[10]*(mobility_matrix[i,t-1]-mmobility_matrix)
      lp21op22[i,t] <- alpha[2]+alpha[6]*(standard_cap[i]-mean(standard_cap[1:N]))+
        alpha[11]*(mobility_matrix[i,t-1]-mmobility_matrix)
      #mobility_matrix is already lagged by 3
      #therefore with t-1 it is lagged by 4
      lp23op22[i,t] <- alpha[3]+alpha[7]*(mobility_matrix[i,t-1]-mmobility_matrix)+
        alpha[8]*(standard_cap[i]-mean(standard_cap[1:N]))+
        alpha[13]*new_variant[t]
      p22[i,t] <- 1/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      p21[i,t] <- exp(lp21op22[i,t])/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      p23[i,t] <- exp(lp23op22[i,t])/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      logit(p33[i,t]) <- alpha[4]+alpha[9]*(standard_cap[i]-mean(standard_cap[1:N]))+
        alpha[12]*(mobility_matrix[i,t-1]-mmobility_matrix)
      
      #transit[i,t,1:3] <- tm[i,t,S[i,t-1],1:3]
      S[i,t] ~ dcat(tm[i,t,S[i,t-1],1:7])
    }
  }
  
  #random effects
  for(i in 1:N){
    b0[i] ~ dnorm(mean=beta0+beta2*(standard_cap[i]-mean(standard_cap[1:N])),prec_b0)
    b1[i] ~ dnorm(mean=beta1+beta3*(standard_cap[i]-mean(standard_cap[1:N])),prec_b1)
  }
  
  #constraints
  constraint_data1 ~ dconstraint(b0[1]+.1 < b1[1] & b0[2]+.1 < b1[2] & b0[3]+.1 < b1[3] & b0[4]+.1 < b1[4] & b0[5]+.1 < b1[5] & 
                                   b0[6]+.1 < b1[6] & b0[7]+.1 < b1[7] & b0[8]+.1 < b1[8] & b0[9]+.1 < b1[9] & b0[10]+.1 < b1[10])
  
  constraint_data2 ~ dconstraint(b0[11]+.1 < b1[11] & b0[12]+.1 < b1[12] & b0[13]+.1 < b1[13] & b0[14]+.1 < b1[14] & b0[15]+.1 < b1[15] & 
                                   b0[16]+.1 < b1[16] & b0[17]+.1 < b1[17] & b0[18]+.1 < b1[18] & b0[19]+.1 < b1[19] & b0[20]+.1 < b1[20])
  
  constraint_data3 ~ dconstraint(b0[21]+.1 < b1[21] & b0[22]+.1 < b1[22] & b0[23]+.1  < b1[23]& b0[24]+.1 < b1[24] & b0[25]+.1 < b1[25] & 
                                   b0[26]+.1 < b1[26] & b0[27]+.1 < b1[27] & b0[28]+.1 < b1[28] & b0[29]+.1 < b1[29] & b0[30]+.1 < b1[30])
  
  constraint_data4 ~ dconstraint(rho0+.05<rho1)
  
  #priors
  beta0~dnorm(0,sd=100)
  beta1~dnorm(0,sd=100)
  beta2~dnorm(0,sd=100)
  beta3~dnorm(0,sd=100)
  rho1 ~ dunif(.1,1)
  rho0 ~ dunif(.1,1)
  r0 ~ dunif(0,10)
  r1 ~ dunif(0,50)
  prec_b0 ~ dgamma(.1,.1)
  prec_b1 ~ dgamma(.1,.1)
  sigma_b0 <- 1/sqrt(prec_b0)
  sigma_b1 <- 1/sqrt(prec_b1)
  alpha[1] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[2] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[3] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[4] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[5] ~ dt(mu=0,tau=1/(.9442955^2),df=1)
  alpha[6] ~ dt(mu=0,tau=1/(.9442955^2),df=1)
  alpha[7] ~ dt(mu=0,tau=1/(.06017889^2),df=1)
  alpha[8] ~ dt(mu=0,tau=1/(.9442955^2),df=1)
  alpha[9] ~ dt(mu=0,tau=1/(.9442955^2),df=1)
  alpha[10] ~ dt(mu=0,tau=1/(.06017889^2),df=1)
  alpha[11] ~ dt(mu=0,tau=1/(.06017889^2),df=1)
  alpha[12] ~ dt(mu=0,tau=1/(.06017889^2),df=1)
  alpha[13] ~ dt(mu=0,tau=1/(2.5^2),df=1)
  
  
})


#Now we generate valid initial values for the hidden Markov chain
tm_S_init <- matrix(nrow=6,ncol=6)
tm_S_init[1,] <- c(0,1,0,0,0,0)
tm_S_init[2,] <- c(0,.8,.2,0,0,0)
tm_S_init[3,] <- c(0,0,0,1,0,0)
tm_S_init[4,] <- c(0,0,0,0,1,0)
tm_S_init[5,] <- c(0,0,0,0,0,1)
tm_S_init[6,] <- c(.2,0,0,0,0,.8)

S_init <- matrix(nrow=N,ncol=mT)

for(i in 1:N){
  S_init[i,1] <- sample(x=c(1,2,3,4,5,6),size = 1,prob=c(1/6,1/6,1/6,1/6,1/6,1/6))
  for(t in 2:mT){
    S_init[i,t] <- sample(x=c(1,2,3,4,5,6),size = 1,prob=tm_S_init[S_init[i,t-1],])
  }
}

sum(is.na(S_init))

S_init <- S_init+1

#now we generate valid random initial values for the parameters
b0_init <- rnorm(n=N,mean=0,sd=.5) 
rho0_init <- runif(n=1,min=.1,max=.7)
hospitalInits <- list("beta0"=rnorm(n=1,mean=0,sd=.5),"beta1"=rnorm(n=1,mean=0,sd=.5),
                      "beta2"=rnorm(n=1,mean=0,sd=.5),"beta3"=rnorm(n=1,mean=0,sd=.5),
                      "b0"=b0_init,"b1"=b0_init+runif(n=N,min=.1,max=2),
                      "rho0"=rho0_init, "rho1"=runif(n=1,rho0_init+.05,1),
                      "r1"=runif(n=1,min=0,max=20), "r0"=runif(n=1,min=0,max=10),
                      "prec_b0"=runif(n=1,min=0,max=20),"prec_b1"=runif(n=1,min=0,max=20),
                      "alpha"=rnorm(n=13,mean=0,sd=c(rep(.5,9),.1,.1,.1,.1)),
                      "S"=S_init)

#now we build and compile the model
hospital_model <- nimbleModel(code = hospitalCode,inits = hospitalInits, constants = hospitalConsts,
                              data = hospitalData)

chospital_model  <- compileNimble(hospital_model)

#make sure no NAs, if there are NAs then most likely need to tighten the inital values
hospital_model$getLogProb()



#the following code will setup and run an iteration of the FFBS sampler in a single area
#this is very useful for testing/debugging the sampler
model <- hospital_model
dq1 <- !hospital_model$isData("S[1, ]")
target <- "S[1, dq1]"

model$S[1, dq1]

#setup
#nnames
nnames <- model$expandNodeNames(target)
times <- as.numeric(gsub(".*\\[.*,(.*)\\].*", "\\1", nnames))
loc <- as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", nnames[[1]]))
calcNodes <- model$getDependencies(target)
numnodes <- length(nnames)
#grab the dependencies of each target node
dependencies <- NULL
start_depend <- rep(NA,numnodes)
start_depend[1] <- 1
end_depend <- rep(NA,numnodes)
index <- 1
for (n in nnames){
  d <- model$getDependencies(n)
  dependencies <- c(dependencies,d)
  end_depend[index] <- length(d)+start_depend[index]-1
  start_depend[index+1] <- end_depend[index]+1
  index <- index+1
}
q <- matrix(nrow=numnodes,ncol=7)
q[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)
#log likelihood for the states
ll <- matrix(nrow=numnodes,ncol=7)
ll[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)

#now run 
#start with ct=1
#in a non coupled filter it is just the initial state distribution
q[1,1:7] <- model$uniff

#now run filter 
for(ct in 2:numnodes){
  
  #ct <- 2
  q_tm1 <- q[ct-1,1:7]
  p <- t(model$tm[loc,ct,1:7,1:7]) %*% asCol(q_tm1[1:7])
  
  #now need to calculate log-likelood
  
  #state 1 
  if(model$y[loc,ct]==0){
    ll[ct,1] <- 0
  }else{
    ll[ct,1] <- -Inf
  }
  
  #state 2 and 3
  llen <- dnbinom(x=model$y[loc,ct],size = model$r0 ,
                  prob = model$r0/(model$r0+model$lambda_en[loc,ct]),log =TRUE)
  
  ll[ct,2] <- llen
  ll[ct,3] <- llen
  
  #states 4 and 5
  llep <- dnbinom(x=model$y[loc,ct],size = model$r1 ,
                  prob = model$r1/(model$r1+model$lambda_ep[loc,ct]),log =TRUE)
  
  ll[ct,4] <- llep
  ll[ct,5] <- llep
  ll[ct,6] <- llep
  ll[ct,7] <- llep
  
  nl <- ll[ct,1:7]+log(p[1:7,1])
  nls <- nl-max(nl)
  q[ct,1:7] <- exp(nls)/sum(exp(nls))
  
}

#now backward sampler
#first from its filtered prob
prev <- model$S[loc,numnodes]
news <- rcat(1,q[numnodes,1:7])
model$S[loc,numnodes] <- news

if(model$S[loc,numnodes]!= prev){
  model$calculate(nodes = dependencies[start_depend[numnodes]:end_depend[numnodes]])
}

for(ict in 2:numnodes){
  ct <- numnodes-ict+1
  prev <- model$S[loc,ct]
  fs <- model$S[loc,ct+1]
  trans <- model$tm[loc,ct+1,1:7,fs]
  lp <- log(trans)+log(q[ct,1:7])
  lp <- lp-max(lp)
  b <- exp(lp)/sum(exp(lp))
  news <- rcat(1,b)
  model$S[loc,ct] <- news
  if(model$S[loc,ct]!= prev){
    model$calculate(nodes = dependencies[start_depend[ct]:end_depend[ct]])
  }
}

#testing
model$S[1, dq1]
model$getLogProb()
model$calculate()

#FFBS sampler code
FFBS <- nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    nnames <- model$expandNodeNames(target)
    times <- as.numeric(gsub(".*\\[.*,(.*)\\].*", "\\1", nnames))
    loc <- as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", nnames[[1]]))
    calcNodes <- model$getDependencies(target)
    numnodes <- length(nnames)
    #grab the dependencies of each target node
    dependencies <- NULL
    start_depend <- rep(NA,numnodes)
    start_depend[1] <- 1
    end_depend <- rep(NA,numnodes)
    index <- 1
    for (n in nnames){
      d <- model$getDependencies(n)
      dependencies <- c(dependencies,d)
      end_depend[index] <- length(d)+start_depend[index]-1
      start_depend[index+1] <- end_depend[index]+1
      index <- index+1
    }
    q <- matrix(nrow=numnodes,ncol=7)
    q[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)
    #log likelihood for the states
    ll <- matrix(nrow=numnodes,ncol=7)
    ll[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)
    
  },
  
  
  run = function() {
    
    #now run 
    #start with ct=1
    #in a non coupled filter it is just the initial state distribution
    q[1,1:7] <<- model$uniff
    
    #now run filter 
    for(ct in 2:numnodes){
      
      #ct <- 2
      q_tm1 <- q[ct-1,1:7]
      p <- t(model$tm[loc,ct,1:7,1:7]) %*% asCol(q_tm1[1:7])
      
      #now need to calculate log-likelood
      
      #state 1 
      if(model$y[loc,ct]==0){
        ll[ct,1] <<- 0
      }else{
        ll[ct,1] <<- -Inf
      }
      
      #state 2 and 3
      llen <- dnbinom(x=model$y[loc,ct],size = model$r0[1] ,
                      prob = model$r0[1]/(model$r0[1]+model$lambda_en[loc,ct]),log =TRUE)
      
      ll[ct,2] <<- llen
      ll[ct,3] <<- llen
      
      #states 4 and 5
      llep <- dnbinom(x=model$y[loc,ct],size = model$r1[1] ,
                      prob = model$r1[1]/(model$r1[1]+model$lambda_ep[loc,ct]),log =TRUE)
      
      ll[ct,4] <<- llep
      ll[ct,5] <<- llep
      ll[ct,6] <<- llep
      ll[ct,7] <<- llep
      
      nl <- ll[ct,1:7]+log(p[1:7,1])
      nls <- nl-max(nl)
      q[ct,1:7] <<- exp(nls)/sum(exp(nls))
      
    }
    
    #now backward sampler
    #first from its filtered prob
    prev <- model$S[loc,numnodes]
    news <- rcat(1,q[numnodes,1:7])
    model$S[loc,numnodes] <<- news
    
    if(model$S[loc,numnodes]!= prev){
      model$calculate(nodes = dependencies[start_depend[numnodes]:end_depend[numnodes]])
    }
    
    for(ict in 2:numnodes){
      ct <- numnodes-ict+1
      prev <- model$S[loc,ct]
      fs <- model$S[loc,ct+1]
      trans <- model$tm[loc,ct+1,1:7,fs]
      lp <- log(trans)+log(q[ct,1:7])
      lp <- lp-max(lp)
      b <- exp(lp)/sum(exp(lp))
      news <- rcat(1,b)
      model$S[loc,ct] <<- news
      if(model$S[loc,ct]!= prev){
        model$calculate(nodes = dependencies[start_depend[ct]:end_depend[ct]])
      }
    }
    
    
    
    copy(from = model, to = mvSaved, row = 1, 
         nodes = calcNodes, logProb = TRUE)
    
  },
  
  methods = list(   reset = function () {}   )
  
)



#now setup the MCMC 
hospital_modelConf <- configureMCMC(hospital_model, print = TRUE)

#add the FFBS samplers
for(loc in 1:N){
  
  dq <- !hospital_model$isData(paste0("S[",loc,", ]"))
  hospital_modelConf$removeSampler(paste0("S[",loc,", dq]"))
  hospital_modelConf$addSampler(target = hospital_model$expandNodeNames(paste0("S[",loc,", dq]")),
                                type = "FFBS")
  
}

print(hospital_modelConf)

hospital_modelConf$addMonitors(c("S","b0","b1","sigma_b0","sigma_b1","r0","r1"))

hospitalMCMC <- buildMCMC(hospital_modelConf)

ChospitalMCMC <- compileNimble(hospitalMCMC, project = hospital_model ,resetFunctions = TRUE)

#each chain will generate intial MCMC values using this function
initsFunction <- function(){
  
  b0_init <- rnorm(n=N,mean=0,sd=.5) 
  rho0_init <- runif(n=1,min=.1,max=.7)
  
  return( list("beta0"=rnorm(n=1,mean=0,sd=.5),"beta1"=rnorm(n=1,mean=0,sd=.5),
               "beta2"=rnorm(n=1,mean=0,sd=.5),"beta3"=rnorm(n=1,mean=0,sd=.5),
               "b0"=b0_init,"b1"=b0_init+runif(n=N,min=.1,max=2),
               "rho0"=rho0_init, "rho1"=runif(n=1,rho0_init+.05,1),
               "r1"=runif(n=1,min=0,max=20), "r0"=runif(n=1,min=0,max=10),
               "prec_b0"=runif(n=1,min=0,max=20),"prec_b1"=runif(n=1,min=0,max=20),
               "alpha"=rnorm(n=13,mean=0,sd=c(rep(.5,9),.1,.1,.1,.1)),
               "S"=S_init))
} 

#custom samplers cannot be run in parrelel
samples <- runMCMC(ChospitalMCMC,  niter =200000,nchains = 3,nburnin=50000
                   ,samplesAsCodaMCMC = TRUE,thin=15,inits = initsFunction)


#check convergence

#should all be less then 1.05
gelman.diag(samples[,c("beta0","beta1","beta2","beta3","r1","sigma_b0","sigma_b1","rho1","alpha[1]","alpha[2]",
                       "b0[1]","alpha[3]","alpha[4]","rho0","r0","alpha[5]","alpha[6]","alpha[7]","alpha[8]","alpha[9]",
                       "alpha[10]","alpha[11]","alpha[12]","alpha[13]",
                       "b1[1]","b0[2]","b1[2]","b0[3]","b1[3]","b0[4]","b1[4]","b0[5]","b1[5]",
                       "b0[6]","b1[6]","b0[7]","b1[7]","b0[8]","b1[8]","b0[9]","b1[9]",
                       "b0[10]","b1[10]","b0[11]","b1[11]","b0[12]","b1[12]","b0[13]","b1[13]",
                       "b0[14]","b1[14]","b0[15]","b1[15]","b0[16]","b1[16]",
                       "b0[17]","b1[17]","b0[18]","b1[18]","b0[19]","b1[19]","b0[20]","b1[20]",
                       "b0[21]","b1[21]","b0[22]","b1[22]","b0[23]","b1[23]",
                       "b0[24]","b1[24]","b0[25]","b1[25]","b0[26]","b1[26]",
                       "b0[27]","b1[27]","b0[28]","b1[28]","b0[29]","b1[29]","b0[30]","b1[30]")])


#should be >1000
min(effectiveSize(samples[,c("beta0","beta1","beta2","beta3","r1","sigma_b0","sigma_b1","rho1","alpha[1]","alpha[2]",
                             "b0[1]","alpha[3]","alpha[4]","rho0","r0","alpha[5]","alpha[6]","alpha[7]","alpha[8]","alpha[9]",
                             "alpha[10]","alpha[11]","alpha[12]","alpha[13]",
                             "b1[1]","b0[2]","b1[2]","b0[3]","b1[3]","b0[4]","b1[4]","b0[5]","b1[5]",
                             "b0[6]","b1[6]","b0[7]","b1[7]","b0[8]","b1[8]","b0[9]","b1[9]",
                             "b0[10]","b1[10]","b0[11]","b1[11]","b0[12]","b1[12]","b0[13]","b1[13]",
                             "b0[14]","b1[14]","b0[15]","b1[15]","b0[16]","b1[16]",
                             "b0[17]","b1[17]","b0[18]","b1[18]","b0[19]","b1[19]","b0[20]","b1[20]",
                             "b0[21]","b1[21]","b0[22]","b1[22]","b0[23]","b1[23]",
                             "b0[24]","b1[24]","b0[25]","b1[25]","b0[26]","b1[26]",
                             "b0[27]","b1[27]","b0[28]","b1[28]","b0[29]","b1[29]","b0[30]","b1[30]")]))

#visually examine the traceplots
plot(samples[,c("beta0","beta1")])
plot(samples[,c("beta2","beta3")])
plot(samples[,c("r0","r1")])
plot(samples[,c("sigma_b0","sigma_b1")])
plot(samples[,c("rho0","rho1")])
plot(samples[,c("alpha[1]","alpha[2]")])
plot(samples[,c("alpha[2]","alpha[3]","alpha[4]")])
plot(samples[,c("alpha[5]","alpha[6]","alpha[7]")])
plot(samples[,c("alpha[8]","alpha[9]")])
plot(samples[,c("alpha[10]","alpha[11]","alpha[12]","alpha[13]")])
plot(samples[,c("b0[1]","b1[1]")])
plot(samples[,c("b0[2]","b1[2]")])
plot(samples[,c("b0[3]","b1[3]")])
plot(samples[,c("b0[4]","b1[4]")])
plot(samples[,c("b0[5]","b1[5]")])
plot(samples[,c("b0[6]","b1[6]")])
plot(samples[,c("b0[7]","b1[7]")])
plot(samples[,c("b0[8]","b1[8]")])
plot(samples[,c("b0[9]","b1[9]")])
plot(samples[,c("b0[10]","b1[10]")])
plot(samples[,c("b0[11]","b1[11]")])
plot(samples[,c("b0[12]","b1[12]")])
plot(samples[,c("b0[13]","b1[13]")])
plot(samples[,c("b0[14]","b1[14]")])
plot(samples[,c("b0[15]","b1[15]")])
plot(samples[,c("b0[16]","b1[16]")])
plot(samples[,c("b0[17]","b1[17]")])
plot(samples[,c("b0[18]","b1[18]")])
plot(samples[,c("b0[19]","b1[19]")])
plot(samples[,c("b0[20]","b1[20]")])
plot(samples[,c("b0[21]","b1[21]")])
plot(samples[,c("b0[22]","b1[22]")])
plot(samples[,c("b0[23]","b1[23]")])
plot(samples[,c("b0[24]","b1[24]")])
plot(samples[,c("b0[25]","b1[25]")])
plot(samples[,c("b0[26]","b1[26]")])
plot(samples[,c("b0[27]","b1[27]")])
plot(samples[,c("b0[28]","b1[28]")])
plot(samples[,c("b0[29]","b1[29]")])
plot(samples[,c("b0[30]","b1[30]")])


#this code will calculate the retrospective fits
samps <- data.frame(rbind(samples[[1]],samples[[2]],samples[[3]]))
fitted <- array(dim = c(N,mT,30000))
fitted_en <- array(dim = c(N,mT,30000))
fitted_ep <- array(dim = c(N,mT,30000))
S_ep <- array(dim = c(N,mT,30000))
S_en <- array(dim = c(N,mT,30000))
S_ab <- array(dim = c(N,mT,30000))

for(i in 1:N){
  print(i)
  for(t in 2:mT){
    
    lambda_epit <- exp(as.numeric(unlist(samps[paste0("b1.",i,".")]))+
                         as.numeric(unlist(samps[paste0("rho1")]))*(lpsi[i,t-1]-mlpsi))
    
    
    lambda_enit <-  exp(as.numeric(unlist(samps[paste0("b0.",i,".")]))+
                          as.numeric(unlist(samps[paste0("rho0")]))*(lpsi[i,t-1]-mlpsi))
    
    
    r0it <- as.numeric(unlist(samps[paste0("r0")]))
    r1it <- as.numeric(unlist(samps[paste0("r1")]))
    
    #I see no other way to do this then to go through the draws in a for loop
    #might be extremely slow though
    sit <- as.numeric(unlist(samps[paste0("S.",i,"..",t,".")]))
    lambdait <- rep(NA,30000)
    indiit <- rep(NA,30000)
    rit <- rep(NA,30000)
    for(j in 1:30000){
      lambdait[j] <- c(1,lambda_enit[j],lambda_enit[j],lambda_epit[j],lambda_epit[j],lambda_epit[j],lambda_epit[j])[sit[j]]
      rit[j] <- c(1,r0it[j],r0it[j],r1it[j],r1it[j],r1it[j],r1it[j])[sit[j]]
      indiit[j] <- c(0,1,1,1,1,1,1)[sit[j]]
    }
    
    pit <- (rit/(lambdait*indiit+rit))-.0000001*(1-indiit)
    
    fitted[i,t,] <- rnbinom(n = 30000,
                            prob=pit,size=rit)
    
    fitted_en[i,t,] <- rnbinom(n = 30000,
                               prob = r0it/(r0it+lambda_enit),size = r0it)
    
    fitted_ep[i,t,] <- rnbinom(n = 30000,
                               prob = r1it/(r1it+lambda_epit),size = r1it)
    
    S_ep[i,t,] <- ifelse(sit==4 | sit==5 | sit==6 | sit==7,1,0) 
    
    S_en[i,t,] <- ifelse(sit==2 | sit==3,1,0) 
    
    S_ab[i,t,] <- ifelse(sit==1,1,0) 
    
  }
}


#will plot the retropspective fit like Figure 4 in the main text
for(i in 1:30){
  
  med_en = apply(fitted_en[i,,],MARGIN = 1,mean)
  upper_en= apply(fitted_en[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.975),na.rm=TRUE))
  lower_en= apply(fitted_en[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.025),na.rm=TRUE))
  
  med_ep = apply(fitted_ep[i,,],MARGIN = 1,mean)
  upper_ep= apply(fitted_ep[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.975),na.rm=TRUE))
  lower_ep= apply(fitted_ep[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.025),na.rm=TRUE))
  
  
  par(mfrow=c(2,1))
  
  plot(1:mT,admissions[i,],type = "p",main=paste0(i),ylim=c(0,max(c(upper_en,upper_ep),na.rm = TRUE)+.25*mean(c(med_en,med_ep),na.rm = TRUE)),
       cex=.5)
  
  
  lines(1:mT,med_ep,col="red")
  lines(1:mT,lower_ep,col="red",lty=2)
  lines(1:mT,upper_ep,col="red",lty=2)
  
  lines(1:mT,med_en,col="blue")
  lines(1:mT,lower_en,col="blue",lty=2)
  lines(1:mT,upper_en,col="blue",lty=2)
  
  #states 
  medsep = apply(S_ep[i,,],MARGIN = 1,mean)
  medsen = apply(S_en[i,,],MARGIN = 1,mean)
  medsab = apply(S_ab[i,,],MARGIN = 1,mean)
  plot(medsep,type="l",col="red",ylim=c(0,1))
  lines(medsen,col="blue")
  lines(medsab,col="green")
  
}







