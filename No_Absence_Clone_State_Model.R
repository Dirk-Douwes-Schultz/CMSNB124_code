#code for the No Absence/Clone State Model from Section 5 of the main text
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


#Nimble constants and data
hospitalConsts <- list(mT=length(admissions[1,]),N=length(admissions[,1]),
                       lpsi=lpsi,mlpsi=mlpsi,
                       standard_cap=standard_cap,
                       mobility_matrix=mobility_matrix,
                       mmobility_matrix=mean(mobility_matrix),
                       new_variant=new_variant,
                       adj=adj,num=num,count=count,
                       W_matrix=W_matrix,
                       near_nei_matrix=near_nei_matrix,
                       distance_weights=distance_weights)


hospitalData <- list(y=admissions,
                     constraint_data1=matrix(rep(1,mT*N),nrow=N,ncol=mT),
                     constraint_data2=1)

#model code, not all variables have the same name as the text
hospitalCode <- nimbleCode({
  
  #likelihood
  indi[1:2] <- c(1,1)
  for(i in 1:N){
    for(t in 2:mT){
      
      #data likelihood
      y[i,t] ~ dnegbin(prob=r[i,t,S[i,t]]/(r[i,t,S[i,t]]+lambda[i,t,S[i,t]]*indi[S[i,t]])-1e-10*(1-indi[S[i,t]]),
                       size=r[i,t,S[i,t]])
      
      lambda[i,t,1:2] <- c(lambda_en[i,t],lambda_ep[i,t])
      r[i,t,1:2] <- c(r0,r1)
      lambda_ep[i,t] <- exp(b1[i]+
                              beta4*(mobility_matrix[i,t-1]-mmobility_matrix)+
                              beta6*new_variant[t]+
                              rho1*(lpsi[i,t-1]))
      lambda_en[i,t] <-  exp(b0[i]+
                               beta5*(mobility_matrix[i,t-1]-mmobility_matrix)+
                               rho0*(lpsi[i,t-1]))
    }
  }
  
  #Markov chain
  uniff[1:2] <- c(1/2,1/2)
  epi_indi[1:2] <- c(0,1)
  for(i in 1:N){
    S[i,1] ~ dcat(uniff[1:2])
    for(t in 2:mT){
      
      #calc neiborhood epidemic sum
      # prev_nei_epi_sum[i,t-1] <-(distance_weights[i,1]*epi_indi[S[near_nei_matrix[i,1],t-1]]+
      #                            distance_weights[i,2]*epi_indi[S[near_nei_matrix[i,2],t-1]]+
      #                            distance_weights[i,3]*epi_indi[S[near_nei_matrix[i,3],t-1]]+
      #                            distance_weights[i,4]*epi_indi[S[near_nei_matrix[i,4],t-1]]+
      #                            distance_weights[i,5]*epi_indi[S[near_nei_matrix[i,5],t-1]])
      
      #more efficient then the above commented code
      Snei[i,(t-1),num[i]] <- epi_indi[S[adj[count[i]],(t-1)]]*W_matrix[i,adj[count[i]]]
      for(j in 2:num[i]){
        Snei[i,(t-1),num[i]-j+1] <- Snei[i,(t-1),num[i]-j+2] + 
          epi_indi[S[adj[count[i]+j-1],(t-1)]]*W_matrix[i,adj[count[i]+j-1]]
      }
      
      #calc transition matrix
      tm[i,t,1,1:2] <- c(1-p23[i,t],p23[i,t])
      tm[i,t,2,1:2] <- c(1-p33[i,t],p33[i,t])
      #mobility_matrix is already lagged by 3
      #therefore with t-1 it is lagged by 4
      logit(p23[i,t]) <- alpha[3]+alpha[7]*(mobility_matrix[i,t-1]-mmobility_matrix)+
        alpha[8]*(standard_cap[i]-mean(standard_cap[1:N]))+
        alpha[13]*Snei[i,(t-1),1]+alpha[17]*new_variant[t]
      logit(p33[i,t]) <- alpha[4]+alpha[9]*(standard_cap[i]-mean(standard_cap[1:N]))+
        alpha[12]*(mobility_matrix[i,t-1]-mmobility_matrix)+
        alpha[14]*Snei[i,(t-1),1]
      
      S[i,t] ~ dcat(tm[i,t,S[i,t-1],1:2])
    }
  }
  
  #random effects
  for(i in 1:N){
    b0[i] ~ dnorm(mean=beta0+beta2*(standard_cap[i]-mean(standard_cap[1:N])),prec_b0)
    b1[i] ~ dnorm(mean=beta1+beta3*(standard_cap[i]-mean(standard_cap[1:N])),prec_b1)
  }
  
  #constraints
  for(i in 1:N){
    for(t in 2:mT){
      constraint_data1[i,t] ~ dconstraint((b0[i]+
                                             beta5*(mobility_matrix[i,t-1]-mmobility_matrix)+.01) < (b1[i]+
                                                                                                       beta4*(mobility_matrix[i,t-1]-mmobility_matrix)+
                                                                                                       beta6*new_variant[t]))
    }
  }
  
  constraint_data2 ~ dconstraint(rho0+.05<rho1)
  
  #mlik for the WAIC, see https://groups.google.com/g/nimble-users/c/Essgt2KsEtc/m/X-BAM9O_AwAJ for the idea
  #this will calculate the partially marginalized density for the WAIC
  #see Section 3.1 of the supplementary material
  for(i in 1:N){
    for(t in 1:mT){
      #the dnorm is just fake it will be replaced with a custom sampler
      mlik[i,t] ~ dnorm(0,1)
    }
  }
  
  #priors
  beta0~dnorm(0,sd=100)
  beta1~dnorm(0,sd=100)
  beta2~dnorm(0,sd=100)
  beta3~dnorm(0,sd=100)
  beta4~dnorm(0,sd=100)
  beta5~dnorm(0,sd=100)
  beta6~dnorm(0,sd=100)
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
  alpha[13] ~ dnorm(0,sd=.36)
  alpha[14] ~ dnorm(0,sd=.36)
  alpha[15] ~ dnorm(0,sd=.36)
  alpha[16] ~ dnorm(0,sd=.36)
  alpha[17] ~ dt(mu=0,tau=1/(2.5^2),df=1)
  
  
})

#Now we produce valid starting values for the Markov chain
tm_S_init <- matrix(nrow=2,ncol=2)
tm_S_init[1,] <- c(.8,.2)
tm_S_init[2,] <- c(.2,.8)

S_init <- matrix(nrow=N,ncol=mT)

for(i in 1:N){
  S_init[i,1] <- 2
  for(t in 2:mT){
    S_init[i,t] <- sample(x=c(1,2),size = 1,prob=tm_S_init[S_init[i,t-1],])
  }
}

sum(is.na(S_init))

#now we produce random valid starting values for the parameters
mmobility_matrix <- mean(mobility_matrix)
start <- 0
while(start==0){
  
  b0_init <- rnorm(n=N,mean=0,sd=.5) 
  rho0_init <- runif(n=1,min=.1,max=.7)
  hospitalInits <- list("beta0"=rnorm(n=1,mean=0,sd=.5),"beta1"=rnorm(n=1,mean=0,sd=.5),
                        "beta2"=rnorm(n=1,mean=0,sd=.5),"beta3"=rnorm(n=1,mean=0,sd=.5),
                        "beta4"=rnorm(n=1,mean=0,sd=.1),"beta5"=rnorm(n=1,mean=0,sd=.1),
                        "beta6"=rnorm(n=1,mean=0,sd=.1),
                        "b0"=b0_init,"b1"=b0_init+runif(n=N,min=.1,max=2),
                        "rho0"=rho0_init, "rho1"=runif(n=1,rho0_init+.05,1),
                        "r1"=runif(n=1,min=0,max=20), "r0"=runif(n=1,min=0,max=10),
                        "prec_b0"=runif(n=1,min=0,max=20),"prec_b1"=runif(n=1,min=0,max=20),
                        "alpha"=rnorm(n=17,mean=0,sd=c(rep(.1,9),.1,.1,.1,.1,.1,.1,.1,.1)),
                        "S"=S_init,
                        "mlik"=matrix(rep(0,mT*N),nrow=N,ncol=mT))
  
  len <- matrix(nrow=N,ncol=mT)
  lep <- matrix(nrow=N,ncol=mT)
  for(i in 1:N){
    for(t in 2:mT){
      len[i,t] <- hospitalInits$b0[i]+
        hospitalInits$beta5*(mobility_matrix[i,t-1]-mmobility_matrix)
      lep[i,t] <- hospitalInits$b1[i]+
        hospitalInits$beta4*(mobility_matrix[i,t-1]-mmobility_matrix)+
        hospitalInits$beta6*new_variant[t]
    }
  }
  ifelse(min(lep-len,na.rm=TRUE)>.01,start <- 1,start <- 0)
  
}

min(lep-len,na.rm=TRUE)

hospital_model <- nimbleModel(code = hospitalCode,inits = hospitalInits, constants = hospitalConsts,
                              data = hospitalData)

chospital_model  <- compileNimble(hospital_model)

#make sure no NAs
hospital_model$getLogProb()
#hospital_model$calculate()


#the below code will test the iFFBS sampler in one area
#it will run a single iteration of the sampler in one area
#this is very useful for debugging the sampler
#loc_test <- 5
loc_test <- 16

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
  d <- d[grepl(paste0(".*\\[(?!",loc,",).*,.*\\].*"),d,perl=TRUE)]
  f_dependencies <- c(f_dependencies,d)
  f_end_depend[index] <- length(d)+f_start_depend[index]-1
  f_start_depend[index+1] <- f_end_depend[index]+1
  index <- index+1
}

q <- matrix(nrow=numnodes,ncol=2)
q[1,1:2] <- c(-.99,-.99)
#log likelihood for the states
ll <- matrix(nrow=numnodes,ncol=2)
ll[1,1:2] <- c(-.99,-.99)
#log forward adjustment
lf <- matrix(nrow=numnodes,ncol=2)
lf[1,1:2] <- c(-.99,-.99)

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
  
  #now need to calcl epi_lf
  model$S[loc,1] <- 2
  epi_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,2] <- epi_lf
  
  #now have to go back, no way around this
  model$S[loc,1] <- original_state
  model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
}else{
  #original epi==1
  
  epi_lf <- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,2] <- epi_lf
  
  #now need to calcl end_lf
  model$S[loc,1] <- 1
  end_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,1] <- end_lf
  
  #now have to go back, no way around this
  model$S[loc,1] <- original_state
  model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
  
  
}

#no ll for the first state
nl <- log(model$uniff[1:2])+lf[1,1:2]
nls <- nl-max(nl)
q[1,1:2] <- exp(nls)/sum(exp(nls))


#now run filter 
for(ct in 2:(numnodes-1)){
  
  #ct <- 2
  q_tm1 <- q[ct-1,1:2]
  p <- t(model$tm[loc,ct,1:2,1:2]) %*% asCol(q_tm1[1:2])
  
  #now need to calculate log-likelood
  
  #state 2 and 3
  llen <- dnbinom(x=model$y[loc,ct],size = model$r0 ,
                  prob = model$r0/(model$r0+model$lambda_en[loc,ct]),log =TRUE)
  
  ll[ct,1] <- llen
  
  #states 4 and 5
  llep <- dnbinom(x=model$y[loc,ct],size = model$r1 ,
                  prob = model$r1/(model$r1+model$lambda_ep[loc,ct]),log =TRUE)
  
  ll[ct,2] <- llep
  
  #calc forward dependencies
  original_state <- model$S[loc,ct]
  original_epi <- model$epi_indi[original_state]
  
  if(original_epi==0){
    
    end_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,1] <- end_lf
    
    #now need to calcl epi_lf
    model$S[loc,ct] <- 2
    epi_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,2] <- epi_lf
    
    #now have to go back, no way around this
    model$S[loc,ct] <- original_state
    model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
  }else{
    #original epi==1
    
    epi_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,2] <- epi_lf
    
    #now need to calcl end_lf
    model$S[loc,ct] <- 1
    end_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,1] <- end_lf
    
    #now have to go back, no way around this
    model$S[loc,ct] <- original_state
    model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    
    
  }
  
  nl <- ll[ct,1:2]+log(p[1:2,1])+lf[ct,1:2]
  nls <- nl-max(nl)
  q[ct,1:2] <- exp(nls)/sum(exp(nls))
  
}

#now final filter, no forward dependencies
ct <- numnodes
q_tm1 <- q[ct-1,1:2]
p <- t(model$tm[loc,ct,1:2,1:2]) %*% asCol(q_tm1[1:2])

#now need to calculate log-likelood

#state 2 and 3
llen <- dnbinom(x=model$y[loc,ct],size = model$r0 ,
                prob = model$r0/(model$r0+model$lambda_en[loc,ct]),log =TRUE)

ll[ct,1] <- llen

#states 4 and 5
llep <- dnbinom(x=model$y[loc,ct],size = model$r1 ,
                prob = model$r1/(model$r1+model$lambda_ep[loc,ct]),log =TRUE)

ll[ct,2] <- llep

nl <- ll[ct,1:2]+log(p[1:2,1])
nls <- nl-max(nl)
q[ct,1:2] <- exp(nls)/sum(exp(nls))


#test
#should be the same as before running the test code
#as nothing has been sampled yet
model$S[loc_test, dq1]
model$getLogProb()
model$calculate()

#now backward sampler
#first from its filtered prob
prev <- model$S[loc,numnodes]
news <- rcat(1,q[numnodes,1:2])
model$S[loc,numnodes] <- news

if(model$S[loc,numnodes]!= prev){
  model$calculate(nodes = dependencies[start_depend[numnodes]:end_depend[numnodes]])
}

for(ict in 2:numnodes){
  ct <- numnodes-ict+1
  prev <- model$S[loc,ct]
  fs <- model$S[loc,ct+1]
  trans <- model$tm[loc,ct+1,1:2,fs]
  lp <- log(trans)+log(q[ct,1:2])
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

#compare q with that produced by the WAIC sampler below. should be the same
print(q)

##########################################################
#the code below tests the WAIC sampler
#this calculates the partially marginalized density for the WAIC
#see Section 3.1 of the supplementary materials
#model and loc_test are already set
target <- paste0("mlik[",loc_test,", ",1:mT,"]")

#setup
nnames0 <- model$expandNodeNames(target)
times <- as.numeric(gsub(".*\\[.*,(.*)\\].*", "\\1", nnames0))
loc <- as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", nnames0[[1]]))
calcNodes <- model$getDependencies(target)
numnodes <- length(nnames0)
nnames <- paste0("S[",loc,", ",1:numnodes,"]")
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
  d <- d[grepl(paste0(".*\\[(?!",loc,",).*,.*\\].*"),d,perl=TRUE)]
  f_dependencies <- c(f_dependencies,d)
  f_end_depend[index] <- length(d)+f_start_depend[index]-1
  f_start_depend[index+1] <- f_end_depend[index]+1
  index <- index+1
}

q <- matrix(nrow=numnodes,ncol=2)
q[1,1:2] <- c(-.99,-.99)
#log likelihood for the states
ll <- matrix(nrow=numnodes,ncol=2)
ll[1,1:2] <- c(-.99,-.99)
#log forward adjustment
lf <- matrix(nrow=numnodes,ncol=2)
lf[1,1:2] <- c(-.99,-.99)

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
  
  #now need to calcl epi_lf
  model$S[loc,1] <- 2
  epi_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,2] <- epi_lf
  
  #now have to go back, no way around this
  model$S[loc,1] <- original_state
  model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
}else{
  #original epi==1
  
  epi_lf <- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,2] <- epi_lf
  
  #now need to calcl end_lf
  model$S[loc,1] <- 1
  end_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
  lf[1,1] <- end_lf
  
  #now have to go back, no way around this
  model$S[loc,1] <- original_state
  model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
  
  
}

#no ll for the first state
nl <- log(model$uniff[1:2])+lf[1,1:2]
nls <- nl-max(nl)
q[1,1:2] <- exp(nls)/sum(exp(nls))

#now run filter 
for(ct in 2:(numnodes-1)){
  
  #ct <- 2
  q_tm1 <- q[ct-1,1:2]
  p <- t(model$tm[loc,ct,1:2,1:2]) %*% asCol(q_tm1[1:2])
  
  #now need to calculate log-likelood
  
  #state 2 and 3
  llen <- dnbinom(x=model$y[loc,ct],size = model$r0 ,
                  prob = model$r0/(model$r0+model$lambda_en[loc,ct]),log =TRUE)
  
  ll[ct,1] <- llen
  
  #states 4 and 5
  llep <- dnbinom(x=model$y[loc,ct],size = model$r1 ,
                  prob = model$r1/(model$r1+model$lambda_ep[loc,ct]),log =TRUE)
  
  ll[ct,2] <- llep
  
  nly <- asRow(exp(ll[ct,1:2])) %*% p
  model$mlik[loc,ct] <- nly
  
  #calc forward dependencies
  original_state <- model$S[loc,ct]
  original_epi <- model$epi_indi[original_state]
  
  if(original_epi==0){
    
    end_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,1] <- end_lf
    
    #now need to calcl epi_lf
    model$S[loc,ct] <- 2
    epi_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,2] <- epi_lf
    
    #now have to go back, no way around this
    model$S[loc,ct] <- original_state
    model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
  }else{
    #original epi==1
    
    epi_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,2] <- epi_lf
    
    #now need to calcl end_lf
    model$S[loc,ct] <- 1
    end_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    lf[ct,1] <- end_lf
    
    #now have to go back, no way around this
    model$S[loc,ct] <- original_state
    model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
    
    
  }
  
  nl <- ll[ct,1:2]+log(p[1:2,1])+lf[ct,1:2]
  nls <- nl-max(nl)
  q[ct,1:2] <- exp(nls)/sum(exp(nls))
  
}

#now final filter, no forward dependencies
ct <- numnodes
q_tm1 <- q[ct-1,1:2]
p <- t(model$tm[loc,ct,1:2,1:2]) %*% asCol(q_tm1[1:2])

#now need to calculate log-likelood


#state 2 and 3
llen <- dnbinom(x=model$y[loc,ct],size = model$r0 ,
                prob = model$r0/(model$r0+model$lambda_en[loc,ct]),log =TRUE)

ll[ct,1] <- llen

#states 4 and 5
llep <- dnbinom(x=model$y[loc,ct],size = model$r1 ,
                prob = model$r1/(model$r1+model$lambda_ep[loc,ct]),log =TRUE)

ll[ct,2] <- llep

nly <- asRow(exp(ll[ct,1:2])) %*% p
model$mlik[loc,ct] <- nly

nl <- ll[ct,1:2]+log(p[1:2,1])
nls <- nl-max(nl)
q[ct,1:2] <- exp(nls)/sum(exp(nls))


model$calculate(nodes=calcNodes)

#test
calcNodes
model$S[loc_test, dq1]
model$mlik[loc_test,]
model$getLogProb()
model$calculate()

#q should be the same as above
print(q)


#mlik for WAIC
#this sampler calculates the partially marginalized density for the WAIC
#see Section 3.1 of the supplementary materials
WAIC_mlik <- nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    #setup
    nnames0 <- model$expandNodeNames(target)
    times <- as.numeric(gsub(".*\\[.*,(.*)\\].*", "\\1", nnames0))
    loc <- as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", nnames0[[1]]))
    calcNodes <- model$getDependencies(target)
    numnodes <- length(nnames0)
    nnames <- paste0("S[",loc,", ",1:numnodes,"]")
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
      d <- d[grepl(paste0(".*\\[(?!",loc,",).*,.*\\].*"),d,perl=TRUE)]
      f_dependencies <- c(f_dependencies,d)
      f_end_depend[index] <- length(d)+f_start_depend[index]-1
      f_start_depend[index+1] <- f_end_depend[index]+1
      index <- index+1
    }
    
    q <- matrix(nrow=numnodes,ncol=2)
    q[1,1:2] <- c(-.99,-.99)
    #log likelihood for the states
    ll <- matrix(nrow=numnodes,ncol=2)
    ll[1,1:2] <- c(-.99,-.99)
    #log forward adjustment
    lf <- matrix(nrow=numnodes,ncol=2)
    lf[1,1:2] <- c(-.99,-.99)
    
    
    
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
      
      #now need to calcl epi_lf
      model$S[loc,1] <<- 2
      epi_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,2] <<- epi_lf
      
      #now have to go back, no way around this
      model$S[loc,1] <<- original_state
      model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
    }else{
      #original epi==1
      
      epi_lf <- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,2] <<- epi_lf
      
      #now need to calcl end_lf
      model$S[loc,1] <<- 1
      end_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,1] <<- end_lf
      
      #now have to go back, no way around this
      model$S[loc,1] <<- original_state
      model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
      
      
    }
    
    #no ll for the first state
    nl <- log(model$uniff[1:2])+lf[1,1:2]
    nls <- nl-max(nl)
    q[1,1:2] <<- exp(nls)/sum(exp(nls))
    
    #now run filter 
    for(ct in 2:(numnodes-1)){
      
      #ct <- 2
      q_tm1 <- q[ct-1,1:2]
      p <- t(model$tm[loc,ct,1:2,1:2]) %*% asCol(q_tm1[1:2])
      
      #now need to calculate log-likelood
      
      
      #state 2 and 3
      llen <- dnbinom(x=model$y[loc,ct],size = model$r0[1] ,
                      prob = model$r0[1]/(model$r0[1]+model$lambda_en[loc,ct]),log =TRUE)
      
      ll[ct,1] <<- llen
      
      #states 4 and 5
      llep <- dnbinom(x=model$y[loc,ct],size = model$r1[1] ,
                      prob = model$r1[1]/(model$r1[1]+model$lambda_ep[loc,ct]),log =TRUE)
      
      ll[ct,2] <<- llep
      
      nly <- asRow(exp(ll[ct,1:2])) %*% p[1:2,1]
      model$mlik[loc,ct] <<- nly[1,1]
      
      #calc forward dependencies
      original_state <- model$S[loc,ct]
      original_epi <- model$epi_indi[original_state]
      
      if(original_epi==0){
        
        end_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,1] <<- end_lf
        
        #now need to calcl epi_lf
        model$S[loc,ct] <<- 2
        epi_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,2] <<- epi_lf
        
        #now have to go back, no way around this
        model$S[loc,ct] <<- original_state
        model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
      }else{
        #original epi==1
        
        epi_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,2] <<- epi_lf
        
        #now need to calcl end_lf
        model$S[loc,ct] <<- 1
        end_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,1] <<- end_lf
        
        #now have to go back, no way around this
        model$S[loc,ct] <<- original_state
        model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        
        
      }
      
      nl <- ll[ct,1:2]+log(p[1:2,1])+lf[ct,1:2]
      nls <- nl-max(nl)
      q[ct,1:2] <<- exp(nls)/sum(exp(nls))
      
    }
    
    #now final filter, no forward dependencies
    q_tm1 <- q[numnodes-1,1:2]
    p <- t(model$tm[loc,numnodes,1:2,1:2]) %*% asCol(q_tm1[1:2])
    
    #now need to calculate log-likelood
    
    
    #state 2 and 3
    llen <- dnbinom(x=model$y[loc,numnodes],size = model$r0[1] ,
                    prob = model$r0[1]/(model$r0[1]+model$lambda_en[loc,numnodes]),log =TRUE)
    
    ll[numnodes,1] <<- llen
    
    #states 4 and 5
    llep <- dnbinom(x=model$y[loc,numnodes],size = model$r1[1] ,
                    prob = model$r1[1]/(model$r1[1]+model$lambda_ep[loc,numnodes]),log =TRUE)
    
    ll[numnodes,2] <<- llep
    
    nly <- asRow(exp(ll[numnodes,1:2])) %*% p[1:2,1]
    model$mlik[loc,numnodes] <<- nly[1,1]
    
    nl <- ll[numnodes,1:2]+log(p[1:2,1])
    nls <- nl-max(nl)
    q[numnodes,1:2] <<- exp(nls)/sum(exp(nls))
    
    
    model$calculate(nodes=calcNodes)
    
    
    
    copy(from = model, to = mvSaved, row = 1, 
         nodes = calcNodes, logProb = TRUE)
    
  },
  
  methods = list(   reset = function () {}   )
  
)


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
      d <- d[grepl(paste0(".*\\[(?!",loc,",).*,.*\\].*"),d,perl=TRUE)]
      f_dependencies <- c(f_dependencies,d)
      f_end_depend[index] <- length(d)+f_start_depend[index]-1
      f_start_depend[index+1] <- f_end_depend[index]+1
      index <- index+1
    }
    
    q <- matrix(nrow=numnodes,ncol=2)
    q[1,1:2] <- c(-.99,-.99)
    #log likelihood for the states
    ll <- matrix(nrow=numnodes,ncol=2)
    ll[1,1:2] <- c(-.99,-.99)
    #log forward adjustment
    lf <- matrix(nrow=numnodes,ncol=2)
    lf[1,1:2] <- c(-.99,-.99)
    
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
      
      #now need to calcl epi_lf
      model$S[loc,1] <<- 2
      epi_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,2] <<- epi_lf
      
      #now have to go back, no way around this
      model$S[loc,1] <<- original_state
      model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
    }else{
      #original epi==1
      
      epi_lf <- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,2] <<- epi_lf
      
      #now need to calcl end_lf
      model$S[loc,1] <<- 1
      end_lf <- model$calculate(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      lf[1,1] <<- end_lf
      
      #now have to go back, no way around this
      model$S[loc,1] <<- original_state
      model$calculate(nodes=f_dependencies[f_start_depend[1]:f_end_depend[1]])
      
      
    }
    
    #no ll for the first state
    nl <- log(model$uniff[1:2])+lf[1,1:2]
    nls <- nl-max(nl)
    q[1,1:2] <<- exp(nls)/sum(exp(nls))
    
    #now run filter 
    for(ct in 2:(numnodes-1)){
      
      #ct <- 2
      q_tm1 <- q[ct-1,1:2]
      p <- t(model$tm[loc,ct,1:2,1:2]) %*% asCol(q_tm1[1:2])
      
      #now need to calculate log-likelood
      
      
      #state 2 and 3
      llen <- dnbinom(x=model$y[loc,ct],size = model$r0[1] ,
                      prob = model$r0[1]/(model$r0[1]+model$lambda_en[loc,ct]),log =TRUE)
      
      ll[ct,1] <<- llen
      
      #states 4 and 5
      llep <- dnbinom(x=model$y[loc,ct],size = model$r1[1] ,
                      prob = model$r1[1]/(model$r1[1]+model$lambda_ep[loc,ct]),log =TRUE)
      
      ll[ct,2] <<- llep
      
      #calc forward dependencies
      original_state <- model$S[loc,ct]
      original_epi <- model$epi_indi[original_state]
      
      if(original_epi==0){
        
        end_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,1] <<- end_lf
        
        #now need to calcl epi_lf
        model$S[loc,ct] <<- 2
        epi_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,2] <<- epi_lf
        
        #now have to go back, no way around this
        model$S[loc,ct] <<- original_state
        model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
      }else{
        #original epi==1
        
        epi_lf <- model$getLogProb(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,2] <<- epi_lf
        
        #now need to calcl end_lf
        model$S[loc,ct] <<- 1
        end_lf <- model$calculate(f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        lf[ct,1] <<- end_lf
        
        #now have to go back, no way around this
        model$S[loc,ct] <<- original_state
        model$calculate(nodes=f_dependencies[f_start_depend[ct]:f_end_depend[ct]])
        
        
      }
      
      nl <- ll[ct,1:2]+log(p[1:2,1])+lf[ct,1:2]
      nls <- nl-max(nl)
      q[ct,1:2] <<- exp(nls)/sum(exp(nls))
      
    }
    
    #now final filter, no forward dependencies
    q_tm1 <- q[numnodes-1,1:2]
    p <- t(model$tm[loc,numnodes,1:2,1:2]) %*% asCol(q_tm1[1:2])
    
    #now need to calculate log-likelood
    
    
    #state 2 and 3
    llen <- dnbinom(x=model$y[loc,numnodes],size = model$r0[1] ,
                    prob = model$r0[1]/(model$r0[1]+model$lambda_en[loc,numnodes]),log =TRUE)
    
    ll[numnodes,1] <<- llen
    
    #states 4 and 5
    llep <- dnbinom(x=model$y[loc,numnodes],size = model$r1[1] ,
                    prob = model$r1[1]/(model$r1[1]+model$lambda_ep[loc,numnodes]),log =TRUE)
    
    ll[numnodes,2] <<- llep
    
    nl <- ll[numnodes,1:2]+log(p[1:2,1])
    nls <- nl-max(nl)
    q[numnodes,1:2] <<- exp(nls)/sum(exp(nls))
    
    #now backward sampler
    #first from its filtered prob
    prev <- model$S[loc,numnodes]
    news <- rcat(1,q[numnodes,1:2])
    model$S[loc,numnodes] <<- news
    
    if(model$S[loc,numnodes]!= prev){
      model$calculate(nodes = dependencies[start_depend[numnodes]:end_depend[numnodes]])
    }
    
    for(ict in 2:numnodes){
      ct <- numnodes-ict+1
      prev <- model$S[loc,ct]
      fs <- model$S[loc,ct+1]
      trans <- model$tm[loc,ct+1,1:2,fs]
      lp <- log(trans)+log(q[ct,1:2])
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
    q <- matrix(nrow=numnodes,ncol=2)
    q[1,1:2] <- c(-.99,-.99)
    #log likelihood for the states
    ll <- matrix(nrow=numnodes,ncol=2)
    ll[1,1:2] <- c(-.99,-.99)
    
  },
  
  
  run = function() {
    
    #now run 
    #start with ct=1
    #in a non coupled filter it is just the initial state distribution
    q[1,1:2] <<- model$uniff
    
    #now run filter 
    for(ct in 2:numnodes){
      
      #ct <- 2
      q_tm1 <- q[ct-1,1:2]
      p <- t(model$tm[loc,ct,1:2,1:2]) %*% asCol(q_tm1[1:2])
      
      #now need to calculate log-likelood
      
      
      #state 2 and 3
      llen <- dnbinom(x=model$y[loc,ct],size = model$r0[1] ,
                      prob = model$r0[1]/(model$r0[1]+model$lambda_en[loc,ct]),log =TRUE)
      
      ll[ct,1] <<- llen
      
      #states 4 and 5
      llep <- dnbinom(x=model$y[loc,ct],size = model$r1[1] ,
                      prob = model$r1[1]/(model$r1[1]+model$lambda_ep[loc,ct]),log =TRUE)
      
      ll[ct,2] <<- llep
      
      nl <- ll[ct,1:2]+log(p[1:2,1])
      nls <- nl-max(nl)
      q[ct,1:2] <<- exp(nls)/sum(exp(nls))
      
    }
    
    #now backward sampler
    #first from its filtered prob
    prev <- model$S[loc,numnodes]
    news <- rcat(1,q[numnodes,1:2])
    model$S[loc,numnodes] <<- news
    
    if(model$S[loc,numnodes]!= prev){
      model$calculate(nodes = dependencies[start_depend[numnodes]:end_depend[numnodes]])
    }
    
    for(ict in 2:numnodes){
      ct <- numnodes-ict+1
      prev <- model$S[loc,ct]
      fs <- model$S[loc,ct+1]
      trans <- model$tm[loc,ct+1,1:2,fs]
      lp <- log(trans)+log(q[ct,1:2])
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
  test_depend <- test_depend[grepl(paste0(".*\\[(?!",loc,",).*,.*\\].*"),test_depend,perl=TRUE)]
  
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

#have to set WAIC samplers
for(loc in 1:N){
  hospital_modelConf$removeSampler(paste0("mlik[",loc,", ",1:mT,"]"))
  hospital_modelConf$addSampler(target=hospital_model$expandNodeNames(paste0("mlik[",loc,", ",1:mT,"]")),
                                type="WAIC_mlik")
}

print(hospital_modelConf)

#make sure to check that WAIC_mlik is at the end of the MCMC order
#addsampler should add them to the end but need to check 
hospital_modelConf$printSamplers(executionOrder = TRUE)

hospital_modelConf$addMonitors(c("S","b0","b1","sigma_b0","sigma_b1","r0","r1",
                                 "mlik"))

hospitalMCMC <- buildMCMC(hospital_modelConf)

ChospitalMCMC <- compileNimble(hospitalMCMC, project = hospital_model ,resetFunctions = TRUE)

#this will generate random valid initial values for the MCMC chains
initsFunction <- function(){
  
  start <- 0
  while(start==0){
    
    b0_init <- rnorm(n=N,mean=0,sd=.5) 
    rho0_init <- runif(n=1,min=.1,max=.7)
    hospitalInits <- list("beta0"=rnorm(n=1,mean=0,sd=.5),"beta1"=rnorm(n=1,mean=0,sd=.5),
                          "beta2"=rnorm(n=1,mean=0,sd=.5),"beta3"=rnorm(n=1,mean=0,sd=.5),
                          "beta4"=rnorm(n=1,mean=0,sd=.1),"beta5"=rnorm(n=1,mean=0,sd=.1),
                          "beta6"=rnorm(n=1,mean=0,sd=.1),
                          "b0"=b0_init,"b1"=b0_init+runif(n=N,min=.1,max=2),
                          "rho0"=rho0_init, "rho1"=runif(n=1,rho0_init+.05,1),
                          "r1"=runif(n=1,min=0,max=20), "r0"=runif(n=1,min=0,max=10),
                          "prec_b0"=runif(n=1,min=0,max=20),"prec_b1"=runif(n=1,min=0,max=20),
                          "alpha"=rnorm(n=17,mean=0,sd=c(rep(.1,9),.1,.1,.1,.1,.1,.1,.1,.1)),
                          "S"=S_init,
                          "mlik"=matrix(rep(0,mT*N),nrow=N,ncol=mT))
    
    len <- matrix(nrow=N,ncol=mT)
    lep <- matrix(nrow=N,ncol=mT)
    for(i in 1:N){
      for(t in 2:mT){
        len[i,t] <- hospitalInits$b0[i]+
          hospitalInits$beta5*(mobility_matrix[i,t-1]-mmobility_matrix)
        lep[i,t] <- hospitalInits$b1[i]+
          hospitalInits$beta4*(mobility_matrix[i,t-1]-mmobility_matrix)+
          hospitalInits$beta6*new_variant[t]
      }
    }
    ifelse(min(lep-len,na.rm=TRUE)>.01,start <- 1,start <- 0)
    
  }
  print(min(lep-len,na.rm=TRUE))
  
  return(hospitalInits)
  
} 

#this will run the MCMC
samples <- runMCMC(ChospitalMCMC,  niter =200000,nchains = 3,nburnin=50000
                   ,samplesAsCodaMCMC = TRUE,thin=35,inits = initsFunction)


##test convergence
#should all be less than 1.05
gelman.diag(samples[,c("beta0","beta1","beta2","beta3","r1","sigma_b0","sigma_b1","rho1",
                       "b0[1]","alpha[3]","alpha[4]","rho0","r0","alpha[7]","alpha[8]","alpha[9]",
                       "alpha[12]","alpha[13]","alpha[14]","alpha[17]",
                       "b1[1]","b0[2]","b1[2]","b0[3]","b1[3]","b0[4]","b1[4]","b0[5]","b1[5]",
                       "b0[6]","b1[6]","b0[7]","b1[7]","b0[8]","b1[8]","b0[9]","b1[9]",
                       "b0[10]","b1[10]","b0[11]","b1[11]","b0[12]","b1[12]","b0[13]","b1[13]",
                       "b0[14]","b1[14]","b0[15]","b1[15]","b0[16]","b1[16]",
                       "b0[17]","b1[17]","b0[18]","b1[18]","b0[19]","b1[19]","b0[20]","b1[20]",
                       "b0[21]","b1[21]","b0[22]","b1[22]","b0[23]","b1[23]",
                       "b0[24]","b1[24]","b0[25]","b1[25]","b0[26]","b1[26]",
                       "b0[27]","b1[27]","b0[28]","b1[28]","b0[29]","b1[29]","b0[30]","b1[30]",
                       "beta4","beta5","beta6")])

#should be greater than 1000
min(effectiveSize(samples[,c("beta0","beta1","beta2","beta3","r1","sigma_b0","sigma_b1","rho1",
                             "b0[1]","alpha[3]","alpha[4]","rho0","r0","alpha[7]","alpha[8]","alpha[9]",
                             "alpha[12]","alpha[13]","alpha[14]","alpha[17]",
                             "b1[1]","b0[2]","b1[2]","b0[3]","b1[3]","b0[4]","b1[4]","b0[5]","b1[5]",
                             "b0[6]","b1[6]","b0[7]","b1[7]","b0[8]","b1[8]","b0[9]","b1[9]",
                             "b0[10]","b1[10]","b0[11]","b1[11]","b0[12]","b1[12]","b0[13]","b1[13]",
                             "b0[14]","b1[14]","b0[15]","b1[15]","b0[16]","b1[16]",
                             "b0[17]","b1[17]","b0[18]","b1[18]","b0[19]","b1[19]","b0[20]","b1[20]",
                             "b0[21]","b1[21]","b0[22]","b1[22]","b0[23]","b1[23]",
                             "b0[24]","b1[24]","b0[25]","b1[25]","b0[26]","b1[26]",
                             "b0[27]","b1[27]","b0[28]","b1[28]","b0[29]","b1[29]","b0[30]","b1[30]",
                             "beta4","beta5","beta6")]))

#can also examine some traceplots
#note not all these are parameters in the model
plot(samples[,c("mlik[1, 90]","mlik[1, 100]")])
plot(samples[,c("alpha[13]","alpha[14]")])
plot(samples[,c("alpha[7]","alpha[12]")])
plot(samples[,c("beta4","beta5","beta6")])
plot(samples[,c("alpha[17]","alpha[2]")])

plot(samples[,c("beta0","beta1")])
plot(samples[,c("beta1","rho1")])
plot(samples[,c("beta2","beta3")])
plot(samples[,c("beta4","beta5")])
plot(samples[,c("rho0","rho1")])
plot(samples[,c("alpha[5]","alpha[10]","alpha[15]")])

plot(samples[,c("b0[1]","b1[1]")])
plot(samples[,c("b0[16]","b1[16]")])

#now can calculate the WAIC from Table 1
samps <- data.frame(rbind(samples[[1]],samples[[2]],samples[[3]]))
lppd <- 0 
pwaic <- 0

for(i in 1:N){
  print(i)
  for(t in 2:mT){
    
    mlikit <- samps[,paste0("mlik.",i,"..",t,".")]
    
    lppd <- lppd + log(mean(mlikit))
    pwaic <- pwaic + var(log(mlikit))
  }
}

waic <- -2*(lppd-pwaic)
#should be 17,545

#fitted values
mmobility_matrix <- mean(mobility_matrix)
fitted <- array(dim = c(N,mT,30000))
fitted_en <- array(dim = c(N,mT,30000))
fitted_ep <- array(dim = c(N,mT,30000))
S_ep <- array(dim = c(N,mT,30000))
S_en <- array(dim = c(N,mT,30000))

for(i in 1:N){
  print(i)
  for(t in 2:mT){
    
    lambda_epit <- exp(as.numeric(unlist(samps[paste0("b1.",i,".")]))+
                         samps[,"beta4"]*(mobility_matrix[i,t-1]-mmobility_matrix)+
                         samps[,"beta6"]*new_variant[t]+
                         as.numeric(unlist(samps[paste0("rho1")]))*(lpsi[i,t-1]))
    
    
    lambda_enit <-  exp(as.numeric(unlist(samps[paste0("b0.",i,".")]))+
                          samps[,"beta5"]*(mobility_matrix[i,t-1]-mmobility_matrix)+
                          as.numeric(unlist(samps[paste0("rho0")]))*(lpsi[i,t-1]))
    
    
    r0it <- as.numeric(unlist(samps[paste0("r0")]))
    r1it <- as.numeric(unlist(samps[paste0("r1")]))
    
    sit <- as.numeric(unlist(samps[paste0("S.",i,"..",t,".")]))
    lambdait <- rep(NA,30000)
    indiit <- rep(NA,30000)
    rit <- rep(NA,30000)
    for(j in 1:30000){
      lambdait[j] <- c(lambda_enit[j],lambda_epit[j])[sit[j]]
      rit[j] <- c(r0it[j],r1it[j])[sit[j]]
      indiit[j] <- c(1,1)[sit[j]]
    }
    
    pit <- (rit/(lambdait*indiit+rit))-.0000001*(1-indiit)
    
    fitted[i,t,] <- rnbinom(n = 30000,
                            prob=pit,size=rit)
    
    fitted_en[i,t,] <- rnbinom(n = 30000,
                               prob = r0it/(r0it+lambda_enit),size = r0it)
    
    fitted_ep[i,t,] <- rnbinom(n = 30000,
                               prob = r1it/(r1it+lambda_epit),size = r1it)
    
    S_ep[i,t,] <- ifelse(sit==2,1,0) 
    
    S_en[i,t,] <- ifelse(sit==1,1,0) 
    
    
  }
}

#plots Figure 5 from the main text for all areas
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
  plot(medsep,type="l",col="red",ylim=c(0,1))
  lines(medsen,col="blue")
  
}
