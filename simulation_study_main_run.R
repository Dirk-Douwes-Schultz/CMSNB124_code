#will run the spatial model on the simulated data set from Section 4 of the main text
#produces and evaluates the retrospective state estimates
#to produce the real-time state estimates you need to fit the model
#up to time T=100,...,120, this code just fits up to T=120
#set working directory to source file location
rm(list = ls())
library(questionr)
library(dplyr)
library(nimble)
library(coda)
library(tidyr)
library(padr)
library(stringr)
library(lsr)

memory.limit(size=1E10)

load("simulations1_spatial.RData")

print(table(Ssim))

plot(ysim[1,])
lines(ysim[1,])
plot(ysim[6,])
lines(ysim[6,])

#beds covariate
standard_cap <- c(772,282,359,380,321,253,531,320,240,376,487,288,421,244,419,104,
                  257,261,564,256,191,544,274,371,373,341,277,293,324,467)/100

#need to set up the neighborhood structure
#5 clusters of 6 areas, so each area has 5 neighbors
cluster <- rep(NA,30)
cluster[1:6] <- 1
cluster[7:12] <- 2
cluster[13:18] <- 3
cluster[19:24] <- 4
cluster[25:30] <- 5

#make a nighborhood matrix
#1 if neighbor 0 if not
#neighbor if in the same cluster
N_matrix <- matrix(rep(0,30*30),nrow=30,ncol=30)
for(i in 1:30){
  N_matrix[i,cluster==cluster[i]] <- 1
}

diag(N_matrix) <- 0
rowSums(N_matrix)

#make a weight matrix
#we will take the weight 1 if in same cluster and 0 else
#so weight matrix and neighborhood matrix is the same
W_matrix <- N_matrix

#the following variables are functions of the neighborhood/weight matrix
#they help with interacting with the neighborhoods
count <- rep(NA,30)
num <- rep(NA,30)
for(i in 1:30){
  num[i] <- sum(N_matrix[i,])
}
count[1] <- 1 
for(i in 1:29){
  count[i+1]  <- count[i]+num[i]
}
adj <- as.numeric(as.carAdjacency(N_matrix)$adj)


#more variables
admissions <- ysim
lpsi <- log(admissions+1)
mlpsi <- mean(lpsi)
mT <- length(admissions[1,])
N <- length(admissions[,1])

#Nimble constants and data
hospitalConsts <- list(mT=length(admissions[1,]),
                       N=length(admissions[,1]),
                       lpsi=lpsi,mlpsi=mlpsi,
                       standard_cap=standard_cap,
                       adj=adj,num=num,count=count,
                       W_matrix=W_matrix)


hospitalData <- list(y=admissions,constraint_data1=rep(1,30),
                     constraint_data2=1)

#spatial model code
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
                              rho1*(lpsi[i,t-1]))
      lambda_en[i,t] <-  exp(beta0+beta2*(standard_cap[i]-mean(standard_cap[1:N]))+
                               rho0*(lpsi[i,t-1]))
    }
  }
  
  #Markov chain
  uniff[1:7] <- c(1/3,1/6,1/6,1/12,1/12,1/12,1/12)
  epi_indi[1:7] <- c(0,0,0,1,1,1,1)
  for(i in 1:N){
    S[i,1] ~ dcat(uniff[1:7])
    for(t in 2:mT){
      
      #calc Neighbor Sum for t=t-1 
      Snei[i,(t-1),num[i]] <- epi_indi[S[adj[count[i]],(t-1)]]*W_matrix[i,adj[count[i]]]
      for(j in 2:num[i]){
        Snei[i,(t-1),num[i]-j+1] <- Snei[i,(t-1),num[i]-j+2] + 
          epi_indi[S[adj[count[i]+j-1],(t-1)]]*W_matrix[i,adj[count[i]+j-1]]
      }
      
      #calc transition matrix
      tm[i,t,1,1:7] <- c(1-p12[i,t],p12[i,t],0,0,0,0,0)
      tm[i,t,2,1:7] <- c(0,0,1,0,0,0,0)
      tm[i,t,3,1:7] <- c(p21[i,t],0,p22[i,t],p23[i,t],0,0,0)
      tm[i,t,4,1:7] <- c(0,0,0,0,1,0,0)
      tm[i,t,5,1:7] <- c(0,0,0,0,0,1,0)
      tm[i,t,6,1:7] <- c(0,0,0,0,0,0,1)
      tm[i,t,7,1:7] <- c(0,1-p33[i,t],0,0,0,0,p33[i,t])
      logit(p12[i,t]) <- alpha[1]
      lp21op22[i,t] <- alpha[2]
      #mobility_matrix is already lagged by 3
      #therefore with t-1 it is lagged by 4
      lp23op22[i,t] <- alpha[3]+alpha[5]*Snei[i,(t-1),1]
      p22[i,t] <- 1/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      p21[i,t] <- exp(lp21op22[i,t])/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      p23[i,t] <- exp(lp23op22[i,t])/(1+exp(lp21op22[i,t])+exp(lp23op22[i,t]))
      logit(p33[i,t]) <- alpha[4]
      
      S[i,t] ~ dcat(tm[i,t,S[i,t-1],1:7])
    }
  }
  
  
  #constraints
  #since there are no space-time covariates
  #we only have to loop from i=1 to N
  #we use .1 as opposed to .01 since we are constraining
  #the average difference in transmission
  for(i in 1:N){
    
    constraint_data1[i] ~ dconstraint((beta0+beta2*(standard_cap[i]-mean(standard_cap[1:N]))+.1) < (beta1+beta3*(standard_cap[i]-mean(standard_cap[1:N]))))
  }
  
  constraint_data2 ~ dconstraint(rho0+.05<rho1)
  
  #priors
  beta0~dnorm(0,sd=100)
  beta1~dnorm(0,sd=100)
  beta2~dnorm(0,sd=100)
  beta3~dnorm(0,sd=100)
  rho1 ~ dunif(.1,1)
  rho0 ~ dunif(.1,1)
  r0 <- r1
  r1 ~ dunif(0,50)
  alpha[1] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[2] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[3] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[4] ~ dt(mu=0,tau=1/(5^2),df=1)
  alpha[5] ~ dnorm(0,sd=.36)
  
})

#Now we produce valid starting values for the Markov chain
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

#now we produce random valid starting values for the parameters
start <- 0
while(start==0){
  
  rho0_init <- runif(n=1,min=.1,max=.7)
  hospitalInits <- list("beta0"=rnorm(n=1,mean=0,sd=.5),"beta1"=rnorm(n=1,mean=0,sd=.5),
                        "beta2"=rnorm(n=1,mean=0,sd=.5),"beta3"=rnorm(n=1,mean=0,sd=.5),
                        "rho0"=rho0_init, "rho1"=runif(n=1,rho0_init+.05,1),
                        "r1"=runif(n=1,min=0,max=20),
                        "alpha"=rnorm(n=5,mean=0,sd=.5),
                        "S"=S_init)
  
  len <- matrix(nrow=N,ncol=mT)
  lep <- matrix(nrow=N,ncol=mT)
  for(i in 1:N){
    for(t in 2:mT){
      len[i,t] <- hospitalInits$beta0+hospitalInits$beta2*(standard_cap[i]-mean(standard_cap[1:N]))
      lep[i,t] <- hospitalInits$beta1+hospitalInits$beta3*(standard_cap[i]-mean(standard_cap[1:N]))
    }
  }
  ifelse(min(lep-len,na.rm=TRUE)>.1,start <- 1,start <- 0)
  
}

min(lep-len,na.rm=TRUE)

hospital_model <- nimbleModel(code = hospitalCode,inits = hospitalInits, constants = hospitalConsts,
                              data = hospitalData)

chospital_model  <- compileNimble(hospital_model)

#make sure no NAs
hospital_model$getLogProb()


#the below code will test the iFFBS sampler in one area
#it will run a single iteration of the sampler in one area
#this is very useful for debugging the sampler
loc_test <- 1

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

hospital_modelConf$addMonitors(c("S","r0","r1"))

hospitalMCMC <- buildMCMC(hospital_modelConf)

ChospitalMCMC <- compileNimble(hospitalMCMC, project = hospital_model ,resetFunctions = TRUE)

#this will generate random valid initial values for the MCMC chains
initsFunction <- function(){
  
  start <- 0
  while(start==0){
    
    rho0_init <- runif(n=1,min=.1,max=.7)
    hospitalInits <- list("beta0"=rnorm(n=1,mean=0,sd=.5),"beta1"=rnorm(n=1,mean=0,sd=.5),
                          "beta2"=rnorm(n=1,mean=0,sd=.5),"beta3"=rnorm(n=1,mean=0,sd=.5),
                          "rho0"=rho0_init, "rho1"=runif(n=1,rho0_init+.05,1),
                          "r1"=runif(n=1,min=0,max=20),
                          "alpha"=rnorm(n=5,mean=0,sd=.5),
                          "S"=S_init)
    
    len <- matrix(nrow=N,ncol=mT)
    lep <- matrix(nrow=N,ncol=mT)
    for(i in 1:N){
      for(t in 2:mT){
        len[i,t] <- hospitalInits$beta0+hospitalInits$beta2*(standard_cap[i]-mean(standard_cap[1:N]))
        lep[i,t] <- hospitalInits$beta1+hospitalInits$beta3*(standard_cap[i]-mean(standard_cap[1:N]))
      }
    }
    ifelse(min(lep-len,na.rm=TRUE)>.1,start <- 1,start <- 0)
    
  }
  print(min(lep-len,na.rm=TRUE))
  
  return(hospitalInits)
} 

#this will run the MCMC
samples <- runMCMC(ChospitalMCMC,  niter =200000,nchains = 3,nburnin=50000
                   ,samplesAsCodaMCMC = TRUE,thin=15,inits = initsFunction)



##check convergence
#should all be less than 1.05
gelman.diag(samples[,c("beta0","beta1",
                       "beta2","beta3",
                       "rho0","rho1","r1",
                       "alpha[1]","alpha[2]","alpha[3]","alpha[4]",
                       "alpha[5]")])

ess <- effectiveSize(samples[,c("beta0","beta1",
                                "beta2","beta3",
                                "rho0","rho1","r1",
                                "alpha[1]","alpha[2]","alpha[3]","alpha[4]",
                                "alpha[5]")])

ess
#should be greater than 1000
min(ess)

#fitted values
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
    
    lambda_epit <- exp(as.numeric(unlist(samps[paste0("beta1")]+
                                           samps[paste0("beta3")]*(standard_cap[i]-mean(standard_cap))))+
                         as.numeric(unlist(samps[paste0("rho1")]))*(lpsi[i,t-1]))
    
    
    lambda_enit <-  exp(as.numeric(unlist(samps[paste0("beta0")]+
                                            samps[paste0("beta2")]*(standard_cap[i]-mean(standard_cap))))+
                          as.numeric(unlist(samps[paste0("rho0")]))*(lpsi[i,t-1]))
    
    
    r0it <- as.numeric(unlist(samps[paste0("r0")]))
    r1it <- as.numeric(unlist(samps[paste0("r1")]))
    
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

library(plotrix)

#this will compare the true and estimated states
#also will plot the fit of y
for(i in 1:30){
  
  med_en = apply(fitted_en[i,,],MARGIN = 1,mean)
  upper_en= apply(fitted_en[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.975),na.rm=TRUE))
  lower_en= apply(fitted_en[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.025),na.rm=TRUE))
  
  med_ep = apply(fitted_ep[i,,],MARGIN = 1,mean)
  upper_ep= apply(fitted_ep[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.975),na.rm=TRUE))
  lower_ep= apply(fitted_ep[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.025),na.rm=TRUE))
  
  
  par(mfrow=c(2,1))
  
  plot(1:mT,ysim[i,],type = "p",main=paste0(i),ylim=c(0,max(c(upper_en,upper_ep),na.rm = TRUE)+.25*mean(c(med_en,med_ep),na.rm = TRUE)),
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
  twoord.plot(1:120,medsep,1:120,Ssim[i,],type=c("l","p"),lcol="red",rcol="black",cex=.75)
  #plot(medsep,type="l",col="red",ylim=c(0,1))
  lines(medsen,col="blue")
  lines(medsab,col="green")
  
}

#now we will calculate the ROC curve based on the true states
#versus the retrospective state estimates
epidemic_actual <- ifelse(Ssim==3,1,0)
medsep_matrix <- matrix(nrow=30,ncol=120)
for(i in 1:30){
  
  medsep_matrix[i,] <- apply(S_ep[i,,],MARGIN = 1,mean)
  
}

library(pROC)
par(mfrow=c(1,1))
roc_score=roc(as.vector(epidemic_actual[,2:120]), as.vector(medsep_matrix[,2:120])) #AUC score
plot(roc_score ,main ="ROC curve Outbreaks")
roc_score
#AUC is .995

coords(roc_score, x=.5, input="threshold", ret=c("threshold",
                                                 "specificity", "sensitivity"))
#so will detect 94.4% of outbreak weeks and throw a false alarm during 2% of the endemic weeks
#note retrospective so real-time will be less impressive (see SM Table 4),
#but it does show the state estimates are reasonable and the model,
#is working well

#can be sensitive to the threshold
#it seems raising the threshold above .5 is not worth it
#specificity is already high
coords(roc_score, x=.6, input="threshold", ret=c("threshold",
                                                 "specificity", "sensitivity"))

#more visualization
#compares actual epidemic state and estimated epidemic state
plot(epidemic_actual[6,])
lines(medsep_matrix[6,],col="red")


#calculate timeliness
#cutoff =.5
timeli <- NULL
for(i in 1:N){
  
  #go period by period
  period <- 1:30
  start <- min(which(epidemic_actual[i,period]==1))+min(period)-1
  ooo <- medsep_matrix[i,start:max(period)]
  timeli <- c(timeli,min(which(ooo>.5)))
  if(min(which(ooo>.5))==6 | min(which(ooo>.5))==5) print(i)
  
  period <- 31:60
  start <- min(which(epidemic_actual[i,period]==1))+min(period)-1
  ooo <- medsep_matrix[i,start:max(period)]
  timeli <- c(timeli,min(which(ooo>.5)))
  if(min(which(ooo>.5))==6 | min(which(ooo>.5))==5) print(i)
  
  period <- 61:90
  start <- min(which(epidemic_actual[i,period]==1))+min(period)-1
  ooo <- medsep_matrix[i,start:max(period)]
  timeli <- c(timeli,min(which(ooo>.5)))
  if(min(which(ooo>.5))==6 | min(which(ooo>.5))==5) print(i)
  
  period <- 91:120
  start <- min(which(epidemic_actual[i,period]==1))+min(period)-1
  ooo <- medsep_matrix[i,start:max(period)]
  timeli <- c(timeli,min(which(ooo>.5)))
  if(min(which(ooo>.5))==6 | min(which(ooo>.5))==5) print(i)
  
}

mean(timeli)
#timeliness = 1.48 weeks 
#this means on average the outbreak is detected 1.48 weeks after it begins
#since we have weekly data the fastest it can be detected is 1 week
#i.e., the lowest the timeliness can be is 1 week



