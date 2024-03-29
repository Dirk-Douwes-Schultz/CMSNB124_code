#code for a single replication of the simulation study from the supplementary materials (SM) Section 2
#uses the strong constraints
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

#bring in data
load("Quebec_hospital_study_data.RData")
mT <- length(admissions[1,])
N <- length(admissions[,1])
lpsi <- log(admissions+1)
mlpsi <- mean(lpsi)


#Nimble constants for the simulation
simConsts <- list(mT=length(admissions[1,]),
                  N=length(admissions[,1]),
                  standard_cap=standard_cap,
                  mobility_matrix=mobility_matrix,
                  mmobility_matrix=mean(mobility_matrix),
                  near_nei_matrix=near_nei_matrix,
                  distance_weights=distance_weights,
                  new_variant=new_variant,
                  actual_y=admissions)

#simulation model code
simCode <- nimbleCode({
  
  #likelihood
  y[1:N,1] <- actual_y[1:N,1]
  indi[1:7] <- c(0,1,1,1,1,1,1)
  for(i in 1:N){
    for(t in 2:mT){
      
      #data likelihood
      y[i,t] ~ dnegbin(prob=r[i,t,S[i,t]]/(r[i,t,S[i,t]]+lambda[i,t,S[i,t]]*indi[S[i,t]])-1e-10*(1-indi[S[i,t]]),
                       size=r[i,t,S[i,t]])
      
      lambda[i,t,1:7] <- c(1,lambda_en[i,t],lambda_en[i,t],lambda_ep[i,t],lambda_ep[i,t],lambda_ep[i,t],lambda_ep[i,t])
      r[i,t,1:7] <- c(1,r0,r0,r1,r1,r1,r1)
      lambda_ep[i,t] <- exp(beta1+beta3*(standard_cap[i]-mean(standard_cap[1:N]))+
                              beta4*(mobility_matrix[i,t-1]-mmobility_matrix)+
                              rho1*log(y[i,t-1]+1))
      lambda_en[i,t] <-  exp(beta0+beta2*(standard_cap[i]-mean(standard_cap[1:N]))+
                               beta5*(mobility_matrix[i,t-1]-mmobility_matrix)+
                               rho0*log(y[i,t-1]+1))
    }
  }
  
  #Markov chain
  uniff[1:7] <- c(1/3,1/6,1/6,1/12,1/12,1/12,1/12)
  epi_indi[1:7] <- c(0,0,0,1,1,1,1)
  for(i in 1:N){
    S[i,1] ~ dcat(uniff[1:7])
    for(t in 2:mT){
      #calc neiborhood epidemic sum
      prev_nei_epi_sum[i,t] <-(distance_weights[i,1]*epi_indi[S[near_nei_matrix[i,1],t-1]]+
                                 distance_weights[i,2]*epi_indi[S[near_nei_matrix[i,2],t-1]]+
                                 distance_weights[i,3]*epi_indi[S[near_nei_matrix[i,3],t-1]]+
                                 distance_weights[i,4]*epi_indi[S[near_nei_matrix[i,4],t-1]]+
                                 distance_weights[i,5]*epi_indi[S[near_nei_matrix[i,5],t-1]])
      #calc transition matrix
      tm[i,t,1,1:7] <- c(1-p12[i,t],p12[i,t],0,0,0,0,0)
      tm[i,t,2,1:7] <- c(0,0,1,0,0,0,0)
      tm[i,t,3,1:7] <- c(p21[i,t],0,p22[i,t],p23[i,t],0,0,0)
      tm[i,t,4,1:7] <- c(0,0,0,0,1,0,0)
      tm[i,t,5,1:7] <- c(0,0,0,0,0,1,0)
      tm[i,t,6,1:7] <- c(0,0,0,0,0,0,1)
      tm[i,t,7,1:7] <- c(0,1-p33[i,t],0,0,0,0,p33[i,t])
      logit(p12[i,t]) <- alpha[1]+alpha[2]*(standard_cap[i]-mean(standard_cap[1:N]))
      lp21op22[i,t] <- alpha[3]+alpha[4]*(standard_cap[i]-mean(standard_cap[1:N]))+
        alpha[12]*(mobility_matrix[i,t-1]-mmobility_matrix)
      #mobility_matrix is already lagged by 3
      #therefore with t-1 it is lagged by 4
      lp23op22[i,t] <- alpha[5]+alpha[6]*(mobility_matrix[i,t-1]-mmobility_matrix)+
        alpha[10]*prev_nei_epi_sum[i,t]+alpha[7]*new_variant[t]
      p22[i,t] <- 1/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      p21[i,t] <- exp(lp21op22[i,t])/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      p23[i,t] <- exp(lp23op22[i,t])/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      logit(p33[i,t]) <- alpha[8]+
        alpha[9]*(mobility_matrix[i,t-1]-mmobility_matrix)+
        alpha[11]*prev_nei_epi_sum[i,t]
      
      S[i,t] ~ dcat(tm[i,t,S[i,t-1],1:7])
    }
  }
  
  #priors
  beta0~dnorm(0,sd=100)
  beta1~dnorm(0,sd=100)
  beta2~dnorm(0,sd=100)
  beta3~dnorm(0,sd=100)
  beta4~dnorm(0,sd=100)
  beta5~dnorm(0,sd=100)
  rho1 ~ dunif(.1,1)
  rho0 ~ dunif(.1,1)
  r1 ~ dunif(0,20)
  r0 <- r1
  for(j in 1:12){
    alpha[j] ~ dnorm(0,sd=5)
  }
  
  
})

#true parameter values
inits_sim <- list(beta0 = 0,beta1 = .78, beta2=.17 ,beta3=.06,beta4=.007,
                  beta5=.003,rho0=.65,rho1=.75,r1=10,
                  alpha=c(-.76,.45,-3.6,-.9,-4.15,.025,2.5,2,.025,1.15,.45,-.035))

#this will generate many warnings because the model is not full initialized yet
#the warnings are not concerning
sim_model <- nimbleModel(code = simCode,inits = inits_sim, constants = simConsts)

csim_model <- compileNimble(sim_model)

#can investigate whether the strong constraints are satisfied
mmobility_matrix <- mean(mobility_matrix)
ltrep <- matrix(nrow=30,ncol=113)
ltren <- matrix(nrow=30,ncol=113)
for(i in 1:N){
  for(t in 2:mT){
    ltrep[i,t] <- sim_model$beta1+sim_model$beta3*(standard_cap[i]-mean(standard_cap[1:N]))+
      sim_model$beta4*(mobility_matrix[i,t-1]-mmobility_matrix)
    ltren[i,t] <- sim_model$beta0+sim_model$beta2*(standard_cap[i]-mean(standard_cap[1:N]))+
      sim_model$beta5*(mobility_matrix[i,t-1]-mmobility_matrix)
  }
}
mean(ltrep-ltren,na.rm=TRUE)
min(ltrep-ltren,na.rm=TRUE)
which((ltrep-ltren)==min(ltrep-ltren,na.rm=TRUE),arr.ind = TRUE)

#check values pre-simulation
sim_model$alpha
sim_model$y
sim_model$y[1:30,1]
sim_model$S
sim_model$beta0
sim_model$beta1
sim_model$beta2
sim_model$beta3
sim_model$beta4
sim_model$beta5
sim_model$rho0
sim_model$rho1
sim_model$r1
sim_model$r0

nodesToSim <- csim_model$getDependencies(c("beta0","beta1" ,"beta2","beta3","beta4",
                                           "beta5","alpha","r1",
                                           "rho0","rho1"),
                                         self = F, downstream = T)
#have to add "S" on there
nodesToSim <- c(nodesToSim,"S","prev_nei_epi_sum")

#will simulate from the model
csim_model$simulate(nodesToSim)

#check values post-simulation
csim_model$alpha
csim_model$y
csim_model$S
csim_model$beta0
csim_model$beta1
csim_model$beta2
csim_model$beta3
csim_model$rho0
csim_model$rho1
csim_model$r1

#can plot some of the simulated time series
par(mfrow=c(2,1))
plot(csim_model$y[1,])
lines(csim_model$y[1,])
plot(csim_model$S[1,])

plot(csim_model$y[16,])
lines(csim_model$y[16,])
plot(csim_model$S[16,])

plot(csim_model$y[6,])
lines(csim_model$y[6,])
plot(csim_model$S[6,])

plot(csim_model$y[5,])
lines(csim_model$y[5,])
plot(csim_model$S[5,])

simy <- csim_model$y
simS <- csim_model$S

#end of sims
#############################
#############################

lpsi <- log(csim_model$y+1)

#this should produce a valid starting chain for S
tm_S_init <- matrix(nrow=6,ncol=6)
tm_S_init[1,] <- c(0,1,0,0,0,0)
tm_S_init[2,] <- c(0,.8,.2,0,0,0)
tm_S_init[3,] <- c(0,0,0,1,0,0)
tm_S_init[4,] <- c(0,0,0,0,1,0)
tm_S_init[5,] <- c(0,0,0,0,0,1)
tm_S_init[6,] <- c(.2,0,0,0,0,.8)

S_init <- matrix(nrow=N,ncol=mT)

for(i in 1:N){
  
  S_init[i,1] <- 2
  for(t in 2:mT){
    S_init[i,t] <- sample(x=c(1,2,3,4,5,6),size = 1,prob=tm_S_init[S_init[i,t-1],])
  }
}

sum(is.na(S_init))

S_init <- S_init+1


#model code for the model fit to the simulated data
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
      lambda_ep[i,t] <- exp(beta1+beta3*(standard_cap[i]-mean(standard_cap[1:N]))+
                              beta4*(mobility_matrix[i,t-1]-mmobility_matrix)+
                              rho1*lpsi[i,t-1])
      lambda_en[i,t] <-  exp(beta0+beta2*(standard_cap[i]-mean(standard_cap[1:N]))+
                               beta5*(mobility_matrix[i,t-1]-mmobility_matrix)+
                               rho0*lpsi[i,t-1])
    }
  }
  
  #Markov chain
  uniff[1:7] <- c(1/3,1/6,1/6,1/12,1/12,1/12,1/12)
  epi_indi[1:7] <- c(0,0,0,1,1,1,1)
  for(i in 1:N){
    S[i,1] ~ dcat(uniff[1:7])
    for(t in 2:mT){
      #calc neiborhood epidemic sum
      prev_nei_epi_sum[i,t] <-(distance_weights[i,1]*epi_indi[S[near_nei_matrix[i,1],t-1]]+
                                 distance_weights[i,2]*epi_indi[S[near_nei_matrix[i,2],t-1]]+
                                 distance_weights[i,3]*epi_indi[S[near_nei_matrix[i,3],t-1]]+
                                 distance_weights[i,4]*epi_indi[S[near_nei_matrix[i,4],t-1]]+
                                 distance_weights[i,5]*epi_indi[S[near_nei_matrix[i,5],t-1]])
      #calc transition matrix
      tm[i,t,1,1:7] <- c(1-p12[i,t],p12[i,t],0,0,0,0,0)
      tm[i,t,2,1:7] <- c(0,0,1,0,0,0,0)
      tm[i,t,3,1:7] <- c(p21[i,t],0,p22[i,t],p23[i,t],0,0,0)
      tm[i,t,4,1:7] <- c(0,0,0,0,1,0,0)
      tm[i,t,5,1:7] <- c(0,0,0,0,0,1,0)
      tm[i,t,6,1:7] <- c(0,0,0,0,0,0,1)
      tm[i,t,7,1:7] <- c(0,1-p33[i,t],0,0,0,0,p33[i,t])
      logit(p12[i,t]) <- alpha[1]+alpha[2]*(standard_cap[i]-mean(standard_cap[1:N]))
      lp21op22[i,t] <- alpha[3]+alpha[4]*(standard_cap[i]-mean(standard_cap[1:N]))+
        alpha[12]*(mobility_matrix[i,t-1]-mmobility_matrix)
      #mobility_matrix is already lagged by 3
      #therefore with t-1 it is lagged by 4
      lp23op22[i,t] <- alpha[5]+alpha[6]*(mobility_matrix[i,t-1]-mmobility_matrix)+
        alpha[10]*prev_nei_epi_sum[i,t]+alpha[7]*new_variant[t]
      p22[i,t] <- 1/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      p21[i,t] <- exp(lp21op22[i,t])/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      p23[i,t] <- exp(lp23op22[i,t])/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      logit(p33[i,t]) <- alpha[8]+
        alpha[9]*(mobility_matrix[i,t-1]-mmobility_matrix)+
        alpha[11]*prev_nei_epi_sum[i,t]
      
      S[i,t] ~ dcat(tm[i,t,S[i,t-1],1:7])
    }
  }
  
  #strong constraints
  constraint_data2 ~ dconstraint((rho0+.05)<rho1)
  for(i in 1:N){
    for(t in 2:mT){
      constraint_data1[i,t] ~ dconstraint((beta0+beta2*(standard_cap[i]-mean(standard_cap[1:N]))+
                                             beta5*(mobility_matrix[i,t-1]-mmobility_matrix)+.01) < (beta1+beta3*(standard_cap[i]-mean(standard_cap[1:N]))+
                                                                                                       beta4*(mobility_matrix[i,t-1]-mmobility_matrix)))
    }
  }
  
  
  #priors
  beta0~dnorm(0,sd=100)
  beta1~dnorm(0,sd=100)
  beta2~dnorm(0,sd=100)
  beta3~dnorm(0,sd=100)
  beta4~dnorm(0,sd=100)
  beta5~dnorm(0,sd=100)
  rho1 ~ dunif(.1,1)
  rho0 ~ dunif(.1,1)
  r1 ~ dunif(0,20)
  r0 <- r1
  alpha[1] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[2] ~ dt(mu=0,tau=1/(.9442955^2),df=1)
  alpha[3] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[4] ~ dt(mu=0,tau=1/(.9442955^2),df=1)
  alpha[5] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[6] ~ dt(mu=0,tau=1/(.06017889^2),df=1)
  alpha[7] ~ dt(mu=0,tau=1/(2.5^2),df=1)
  alpha[8] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[9] ~ dt(mu=0,tau=1/(.06017889^2),df=1)
  alpha[10] ~ dt(mu=0,tau=1/(1.25^2),df=1)
  alpha[11] ~ dt(mu=0,tau=1/(1.25^2),df=1)
  alpha[12] ~ dt(mu=0,tau=1/(1.25^2),df=1)
  
  
})

mmobility_matrix <- mean(mobility_matrix)

#model constants and data
hospitalConsts <- list(mT=length(admissions[1,]),
                       N=length(admissions[,1]),lpsi=log(simy+1),
                       standard_cap=standard_cap,
                       mobility_matrix=mobility_matrix,
                       mmobility_matrix=mean(mobility_matrix),
                       near_nei_matrix=near_nei_matrix,
                       distance_weights=distance_weights,
                       new_variant=new_variant)


#note we pass the simulated data as the outcome variable y
hospitalData <- list(y=simy,constraint_data1=matrix(rep(1,mT*N),nrow=N,ncol=mT),
                     constraint_data2=1)

#produces random valid starting values for the parameters
start <- 0
while(start==0){
  b0_init <- rnorm(n=1,mean=0,sd=.5) 
  rho0_init <- runif(n=1,min=.1,max=.7)
  hospitalInits <- list("beta0"=b0_init,"beta1"=b0_init+runif(n=1,min=.1,max=2),
                        "beta2"=rnorm(n=1,mean=0,sd=.5),"beta3"=rnorm(n=1,mean=0,sd=.5),
                        "beta4"=rnorm(n=1,mean=0,sd=.1),"beta5"=rnorm(n=1,mean=0,sd=.1),
                        "rho0"=rho0_init, "rho1"=runif(n=1,rho0_init+.05,1),
                        "r1"=runif(n=1,min=0,max=20),
                        "alpha"=rnorm(n=12,mean=0,sd=.1),
                        "S"=S_init)
  len <- matrix(nrow=N,ncol=mT)
  lep <- matrix(nrow=N,ncol=mT)
  for(i in 1:N){
    for(t in 2:mT){
      len[i,t] <- hospitalInits$beta0+hospitalInits$beta2*(standard_cap[i]-mean(standard_cap[1:N]))+
        hospitalInits$beta5*(mobility_matrix[i,t-1]-mmobility_matrix)
      lep[i,t] <- hospitalInits$beta1+hospitalInits$beta3*(standard_cap[i]-mean(standard_cap[1:N]))+
        hospitalInits$beta4*(mobility_matrix[i,t-1]-mmobility_matrix)
    }
  }
  ifelse(min(lep-len,na.rm=TRUE)>.01,start <- 1,start <- 0)
}

min(lep-len,na.rm = TRUE)


hospital_model <- nimbleModel(code = hospitalCode,inits = hospitalInits, constants = hospitalConsts,
                              data = hospitalData)

chospital_model  <- compileNimble(hospital_model)

#make sure no NAs
hospital_model$getLogProb()


#the below code will test the iFFBS sampler in one area
#it will run a single iteration of the sampler in one area
#this is very useful for debugging the sampler
loc_test <- 5

#arguments
model <- hospital_model
dq1 <- !hospital_model$isData(paste0("S[",loc_test,", ]"))
target <- paste0("S[",loc_test,", dq1]")

model$S[loc_test, dq1]

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
#need forward dependencies, this may have to be modified if you change model specification
#note this does not work for last time as last time has no forward dependencies
f_dependencies <- NULL
f_start_depend <- rep(NA,numnodes)
f_start_depend[1] <- 1
f_end_depend <- rep(NA,numnodes)
index <- 1
for(n in nnames){
  d <- model$getDependencies(n)
  d <- d[grepl(paste0(".*\\[(?!",loc,").*,.*\\].*"),d,perl=TRUE)]
  f_dependencies <- c(f_dependencies,d)
  f_end_depend[index] <- length(d)+f_start_depend[index]-1
  f_start_depend[index+1] <- f_end_depend[index]+1
  index <- index+1
}

q <- matrix(nrow=numnodes,ncol=7)
q[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)
#log likelihood for the states
ll <- matrix(nrow=numnodes,ncol=7)
ll[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)
#log forward adjustment
lf <- matrix(nrow=numnodes,ncol=7)
lf[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)

#now run 
#start with ct=1
#in coupled filter have to add the lf to the log of the initial state distribution
#the forward dependencies should only depend on if the main chain state is epidemic

#calc forward dependencies
original_state <- model$S[loc,1]
original_epi <- model$epi_indi[original_state]

if(original_epi==0){
  
  end_lf <- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,1] <- end_lf
  lf[1,2] <- end_lf
  lf[1,3] <- end_lf
  
  #now need to calcl epi_lf
  model$S[loc,1] <- 4
  epi_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,4] <- epi_lf
  lf[1,5] <- epi_lf
  lf[1,6] <- epi_lf
  lf[1,7] <- epi_lf
  
  #now have to go back, no way around this
  model$S[loc,1] <- original_state
  model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
}else{
  #original epi==1
  
  epi_lf <- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,4] <- epi_lf
  lf[1,5] <- epi_lf
  lf[1,6] <- epi_lf
  lf[1,7] <- epi_lf
  
  #now need to calcl end_lf
  model$S[loc,1] <- 2
  end_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,1] <- end_lf
  lf[1,2] <- end_lf
  lf[1,3] <- end_lf
  
  #now have to go back, no way around this
  model$S[loc,1] <- original_state
  model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
  
  
}

#no ll for the first state
nl <- log(model$uniff[1:7])+lf[1,1:7]
nls <- nl-max(nl)
q[1,1:7] <- exp(nls)/sum(exp(nls))


#now run filter 
for(ct in 2:(numnodes-1)){
  
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
  
  #calc forward dependencies
  original_state <- model$S[loc,ct]
  original_epi <- model$epi_indi[original_state]
  
  if(original_epi==0){
    
    end_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,1] <- end_lf
    lf[ct,2] <- end_lf
    lf[ct,3] <- end_lf
    
    #now need to calcl epi_lf
    model$S[loc,ct] <- 4
    epi_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,4] <- epi_lf
    lf[ct,5] <- epi_lf
    lf[ct,6] <- epi_lf
    lf[ct,7] <- epi_lf
    
    #now have to go back, no way around this
    model$S[loc,ct] <- original_state
    model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
  }else{
    #original epi==1
    
    epi_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,4] <- epi_lf
    lf[ct,5] <- epi_lf
    lf[ct,6] <- epi_lf
    lf[ct,7] <- epi_lf
    
    #now need to calcl end_lf
    model$S[loc,ct] <- 2
    end_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,1] <- end_lf
    lf[ct,2] <- end_lf
    lf[ct,3] <- end_lf
    
    #now have to go back, no way around this
    model$S[loc,ct] <- original_state
    model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    
    
  }
  
  nl <- ll[ct,1:7]+log(p[1:7,1])+lf[ct,1:7]
  nls <- nl-max(nl)
  q[ct,1:7] <- exp(nls)/sum(exp(nls))
  
}

#now final filter, no forward dependencies
ct <- numnodes
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


#test
#should be the same as before running the test code
#as nothing has been sampled yet
model$S[loc_test, dq1]
model$getLogProb()
model$calculate()

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

#now these should be different from before as new states have been sampled
#getLogProb and calculate should return the same values
model$S[loc_test, dq1]
model$getLogProb()
model$calculate()
#so worked properly

#now time to write samplers

#iFFBS
iFFBS <- nimbleFunction(
  
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
    #need forward dependencies, this may have to be modified if you change model specification
    #note this does not work for last time as last time has no forward dependencies
    f_dependencies <- NULL
    f_start_depend <- rep(NA,numnodes)
    f_start_depend[1] <- 1
    f_end_depend <- rep(NA,numnodes)
    index <- 1
    for(n in nnames){
      d <- model$getDependencies(n)
      d <- d[grepl(paste0(".*\\[(?!",loc,").*,.*\\].*"),d,perl=TRUE)]
      f_dependencies <- c(f_dependencies,d)
      f_end_depend[index] <- length(d)+f_start_depend[index]-1
      f_start_depend[index+1] <- f_end_depend[index]+1
      index <- index+1
    }
    
    q <- matrix(nrow=numnodes,ncol=7)
    q[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)
    #log likelihood for the states
    ll <- matrix(nrow=numnodes,ncol=7)
    ll[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)
    #log forward adjustment
    lf <- matrix(nrow=numnodes,ncol=7)
    lf[1,1:7] <- c(-.99,-.99,-.99,-.99,-.99,-.99,-.99)
    
  },
  
  
  run = function() {
    
    #now run 
    #start with ct=1
    #in coupled filter have to add the lf to the log of the initial state distribution
    #the forward dependencies should only depend on if the main chain state is epidemic
    
    #calc forward dependencies
    original_state <- model$S[loc,1]
    original_epi <- model$epi_indi[original_state]
    
    if(original_epi==0){
      
      end_lf <- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,1] <<- end_lf
      lf[1,2] <<- end_lf
      lf[1,3] <<- end_lf
      
      #now need to calcl epi_lf
      model$S[loc,1] <<- 4
      epi_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,4] <<- epi_lf
      lf[1,5] <<- epi_lf
      lf[1,6] <<- epi_lf
      lf[1,7] <<- epi_lf
      
      #now have to go back, no way around this
      model$S[loc,1] <<- original_state
      model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
    }else{
      #original epi==1
      
      epi_lf <- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,4] <<- epi_lf
      lf[1,5] <<- epi_lf
      lf[1,6] <<- epi_lf
      lf[1,7] <<- epi_lf
      
      #now need to calcl end_lf
      model$S[loc,1] <<- 2
      end_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,1] <<- end_lf
      lf[1,2] <<- end_lf
      lf[1,3] <<- end_lf
      
      #now have to go back, no way around this
      model$S[loc,1] <<- original_state
      model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
      
      
    }
    
    #no ll for the first state
    nl <- log(model$uniff[1:7])+lf[1,1:7]
    nls <- nl-max(nl)
    q[1,1:7] <<- exp(nls)/sum(exp(nls))
    
    #now run filter 
    for(ct in 2:(numnodes-1)){
      
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
      
      #calc forward dependencies
      original_state <- model$S[loc,ct]
      original_epi <- model$epi_indi[original_state]
      
      if(original_epi==0){
        
        end_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,1] <<- end_lf
        lf[ct,2] <<- end_lf
        lf[ct,3] <<- end_lf
        
        #now need to calcl epi_lf
        model$S[loc,ct] <<- 4
        epi_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,4] <<- epi_lf
        lf[ct,5] <<- epi_lf
        lf[ct,6] <<- epi_lf
        lf[ct,7] <<- epi_lf
        
        #now have to go back, no way around this
        model$S[loc,ct] <<- original_state
        model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
      }else{
        #original epi==1
        
        epi_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,4] <<- epi_lf
        lf[ct,5] <<- epi_lf
        lf[ct,6] <<- epi_lf
        lf[ct,7] <<- epi_lf
        
        #now need to calcl end_lf
        model$S[loc,ct] <<- 2
        end_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,1] <<- end_lf
        lf[ct,2] <<- end_lf
        lf[ct,3] <<- end_lf
        
        #now have to go back, no way around this
        model$S[loc,ct] <<- original_state
        model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        
        
      }
      
      nl <- ll[ct,1:7]+log(p[1:7,1])+lf[ct,1:7]
      nls <- nl-max(nl)
      q[ct,1:7] <<- exp(nls)/sum(exp(nls))
      
    }
    
    #now final filter, no forward dependencies
    q_tm1 <- q[numnodes-1,1:7]
    p <- t(model$tm[loc,numnodes,1:7,1:7]) %*% asCol(q_tm1[1:7])
    
    #now need to calculate log-likelood
    
    #state 1 
    if(model$y[loc,numnodes]==0){
      ll[numnodes,1] <<- 0
    }else{
      ll[numnodes,1] <<- -Inf
    }
    
    #state 2 and 3
    llen <- dnbinom(x=model$y[loc,numnodes],size = model$r0[1] ,
                    prob = model$r0[1]/(model$r0[1]+model$lambda_en[loc,numnodes]),log =TRUE)
    
    ll[numnodes,2] <<- llen
    ll[numnodes,3] <<- llen
    
    #states 4 and 5
    llep <- dnbinom(x=model$y[loc,numnodes],size = model$r1[1] ,
                    prob = model$r1[1]/(model$r1[1]+model$lambda_ep[loc,numnodes]),log =TRUE)
    
    ll[numnodes,4] <<- llep
    ll[numnodes,5] <<- llep
    ll[numnodes,6] <<- llep
    ll[numnodes,7] <<- llep
    
    nl <- ll[numnodes,1:7]+log(p[1:7,1])
    nls <- nl-max(nl)
    q[numnodes,1:7] <<- exp(nls)/sum(exp(nls))
    
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



#FFBS
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

#now run the model
hospital_modelConf <- configureMCMC(hospital_model, print = TRUE)

#have to add the iFFBS and FFBS samplers
for(loc in 1:N){
  
  #so if the chain for loc has no chains depend on it assign it the FFBS
  # if the chain has other chains depend on it assign the iFFBS
  test_n <- paste0("S[",loc,", 1]")
  test_depend <- model$getDependencies(test_n)
  test_depend <- test_depend[grepl(paste0(".*\\[(?!",loc,").*,.*\\].*"),test_depend,perl=TRUE)]
  
  dq <- !hospital_model$isData(paste0("S[",loc,", ]"))
  hospital_modelConf$removeSampler(paste0("S[",loc,", dq]"))
  
  if(length(test_depend)==0){
    
    hospital_modelConf$addSampler(target = hospital_model$expandNodeNames(paste0("S[",loc,", dq]")),
                                  type = "FFBS")
  }else{
    
    hospital_modelConf$addSampler(target = hospital_model$expandNodeNames(paste0("S[",loc,", dq]")),
                                  type = "iFFBS")
  }
}

print(hospital_modelConf)

hospital_modelConf$addMonitors(c("S","r0","r1"))

hospitalMCMC <- buildMCMC(hospital_modelConf)

ChospitalMCMC <- compileNimble(hospitalMCMC, project = hospital_model ,resetFunctions = TRUE)

#this will generate random valid initial values for the MCMC chains
initsFunction <- function(){
  
  start <- 0
  while(start==0){
    b0_init <- rnorm(n=1,mean=0,sd=.5) 
    rho0_init <- runif(n=1,min=.1,max=.7)
    hospitalInits <- list("beta0"=b0_init,"beta1"=b0_init+runif(n=1,min=.1,max=2),
                          "beta2"=rnorm(n=1,mean=0,sd=.5),"beta3"=rnorm(n=1,mean=0,sd=.5),
                          "beta4"=rnorm(n=1,mean=0,sd=.1),"beta5"=rnorm(n=1,mean=0,sd=.1),
                          "rho0"=rho0_init, "rho1"=runif(n=1,rho0_init+.05,1),
                          "r1"=runif(n=1,min=0,max=20),
                          "alpha"=rnorm(n=12,mean=0,sd=.1),
                          "S"=S_init)
    len <- matrix(nrow=N,ncol=mT)
    lep <- matrix(nrow=N,ncol=mT)
    for(i in 1:N){
      for(t in 2:mT){
        len[i,t] <- hospitalInits$beta0+hospitalInits$beta2*(standard_cap[i]-mean(standard_cap[1:N]))+
          hospitalInits$beta5*(mobility_matrix[i,t-1]-mmobility_matrix)
        lep[i,t] <- hospitalInits$beta1+hospitalInits$beta3*(standard_cap[i]-mean(standard_cap[1:N]))+
          hospitalInits$beta4*(mobility_matrix[i,t-1]-mmobility_matrix)
      }
    }
    ifelse(min(lep-len,na.rm=TRUE)>.01,start <- 1,start <- 0)
  }
  print(min(lep-len,na.rm = TRUE))
  return(hospitalInits)
} 

#this will run the MCMC
samples <- runMCMC(ChospitalMCMC,  niter =200000,nchains = 3,nburnin=50000
                   ,samplesAsCodaMCMC = TRUE,thin=15,inits = initsFunction)


#need to make a summary of the parameters
#note gr=Gelman-Rubin, needs to be less than 1.05 for convergence
#ess=effective sample size, needs to be greater than 1000
#recall the true values from line 116
sim_sum <- data.frame(param=character(0),mean=numeric(0),median=numeric(0),
                      l=numeric(0),u=numeric(0),gr=numeric(0),ess=numeric(0))

params <- c("beta0","beta1","beta2","beta3","beta4","beta5","r1","rho0","rho1","alpha[1]","alpha[2]",
            "alpha[3]","alpha[4]","alpha[5]","alpha[6]","alpha[7]","alpha[8]","alpha[9]",
            "alpha[10]","alpha[11]","alpha[12]")


for(p in params){
  sm <- summary(samples[,c(p)])
  nr <-data.frame(param=p,mean=as.numeric(sm$statistics[1]),
                  median=sm$quantiles[[3]],l=sm$quantiles[[1]],u=sm$quantiles[[5]],
                  gr=gelman.diag(samples[,c(p)])[[1]][[2]],ess=effectiveSize(samples[,c(p)])[[1]])
  sim_sum <- rbind(sim_sum,nr)
}

#can visually examine the traceplots
plot(samples[,c("beta0","beta1")])
plot(samples[,c("beta2","beta3")])
plot(samples[,c("beta4","beta5")])
plot(samples[,c("rho0","rho1","r1","r0")])
plot(samples[,c("alpha[1]","alpha[2]","alpha[3]")])
plot(samples[,c("alpha[4]","alpha[5]","alpha[6]")])
plot(samples[,c("alpha[7]","alpha[8]","alpha[9]")])
plot(samples[,c("alpha[10]","alpha[11]","alpha[12]")])


