#code for the Endemic-epidemic Model from Section 5 of the main text
#fit with a maximum lag of 2, which was the best fitting model
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
#we cleaned the data slightly differently for the EE model
#all data was moved up one week as the model has order 2 autoregression
load("Quebec_hospital_study_data_EE.RData")
mT <- length(admissions[1,])
N <- length(admissions[,1])
lpsi <- log(admissions+1)
mlpsi <- mean(lpsi)
psi <- admissions
mpsi <- mean(admissions)


#Nimble constants and data
hospitalConsts <- list(mT=length(admissions[1,]),
                       N=length(admissions[,1]),
                       lpsi=lpsi,mlpsi=mlpsi,
                       standard_cap=standard_cap,
                       mobility_matrix=mobility_matrix,
                       mmobility_matrix=mean(mobility_matrix),
                       near_nei_matrix=near_nei_matrix,
                       distance_weights=distance_weights,
                       new_variant=new_variant,psi=psi,mpsi=mpsi,
                       adj=adj,count=count,num=num,
                       D_matrix=D_matrix)


hospitalData <- list(y=admissions)

#model code, not all variables are named the same as the main text
hospitalCode <- nimbleCode({
  
  #likelihood
  for(i in 1:N){
    
    #for t=2 need to do it manually
    y_nei_sum[i,(2-1),num[i]] <- psi[adj[count[i]],(2-1)]*W_norm[i,adj[count[i]]]
    for(j in 2:num[i]){
      y_nei_sum[i,(2-1),num[i]-j+1] <- y_nei_sum[i,(2-1),num[i]-j+2] + 
        psi[adj[count[i]+j-1],(2-1)]*W_norm[i,adj[count[i]+j-1]]
    }
    
    for(t in 3:mT){
      
      #calculate the neighborhood sum
      #calc Neighbor Sum for t=t-1 
      y_nei_sum[i,(t-1),num[i]] <- psi[adj[count[i]],(t-1)]*W_norm[i,adj[count[i]]]
      for(j in 2:num[i]){
        y_nei_sum[i,(t-1),num[i]-j+1] <- y_nei_sum[i,(t-1),num[i]-j+2] + 
          psi[adj[count[i]+j-1],(t-1)]*W_norm[i,adj[count[i]+j-1]]
      }
      
      
      y[i,t] ~ dnegbin(prob=r/(r+lambda[i,t]),
                       size=r)
      lambda[i,t] <- exp(b_AR[i]+
                           beta2_AR*(mobility_matrix[i,t-1]-mmobility_matrix)+
                           beta3_AR*new_variant[t])*(y_nei_sum[i,t-1,1]*weight[1]+
                                                       y_nei_sum[i,t-2,1]*weight[2])+
        exp(b_EN[i]+
              beta2_EN*(mobility_matrix[i,t-1]-mmobility_matrix)+
              beta3_EN*new_variant[t])
    }
  }
  
  #spatial weights
  for(i in 1:N){
    for(j in 1:N){
      W[i,j] <- exp(-D_matrix[i,j]*rho)
    }
  }
  #normalize
  for(i in 1:N){
    W_norm[1:N,i] <- W[1:N,i]/sum(W[1:N,i])
  }
  
  
  #random effects
  for(i in 1:N){
    b_AR[i] ~ dnorm(beta0_AR+beta1_AR*(standard_cap[i]-mean(standard_cap[1:N])),prec_AR)
    b_EN[i] ~ dnorm(beta0_EN+beta1_EN*(standard_cap[i]-mean(standard_cap[1:N])),prec_EN) 
  }
  
  #priors
  r~ dunif(0,50)
  rho~dunif(0,1000)
  beta0_AR ~ dnorm(0,100)
  beta0_EN ~ dnorm(0,100)
  prec_AR ~ dgamma(.1,.1)
  prec_EN ~ dgamma(.1,.1)
  beta1_AR ~ dnorm(0,100)
  beta2_AR ~ dnorm(0,100)
  beta3_AR ~ dnorm(0,100)
  beta1_EN ~ dnorm(0,100)
  beta2_EN ~ dnorm(0,100)
  beta3_EN ~ dnorm(0,100)
  
  #lag weights
  for(d in 1:2){
    #geometric weights
    w[d] <- pow(1-kappa,d-1)*kappa
    weight[d] <- w[d]/sum(w[1:2])
    
    #prior distribution of the weights
    w_prior[d] <- pow(1-kappa_prior,d-1)*kappa_prior
    weight_prior[d] <- w_prior[d]/sum(w_prior[1:2])
  }
  kappa ~ dunif(0,1)
  kappa_prior ~ dunif(0,1)
  
  
})

hospitalInits <- list("r"=runif(n=1,min=0,max=20),
                      "beta0_AR"=rnorm(n=1,mean=0,sd=1),
                      "beta0_EN"=rnorm(n=1,mean=0,sd=1),
                      "beta1_AR"=rnorm(n=1,mean=0,sd=1),
                      "beta2_AR"=rnorm(n=1,mean=0,sd=1),
                      "beta3_AR"=rnorm(n=1,mean=0,sd=1),
                      "beta1_EN"=rnorm(n=1,mean=0,sd=1),
                      "beta2_EN"=rnorm(n=1,mean=0,sd=1),
                      "beta3_EN"=rnorm(n=1,mean=0,sd=1),
                      "prec_AR"=runif(n=1,0,10),
                      "prec_EN"=runif(n=1,0,10),
                      "rho"=1,
                      "kappa"=runif(n=1,0,1))

hospital_model <- nimbleModel(code = hospitalCode,inits = hospitalInits, constants = hospitalConsts,
                              data = hospitalData)

chospital_model  <- compileNimble(hospital_model)

hospital_modelConf <- configureMCMC(hospital_model, print = TRUE)

hospital_modelConf$addMonitors(c("b_AR","b_EN","W_norm","weight","weight_prior"))

hospitalMCMC <- buildMCMC(hospital_modelConf)

ChospitalMCMC <- compileNimble(hospitalMCMC, project = hospital_model ,resetFunctions = TRUE)

#this will generate random initial values for the MCMC chains
initsFunction <- function(){
  
  
  return( list("r"=runif(n=1,min=0,max=20),
               "beta0_AR"=rnorm(n=1,mean=0,sd=1),
               "beta0_EN"=rnorm(n=1,mean=0,sd=1),
               "beta1_AR"=rnorm(n=1,mean=0,sd=1),
               "beta2_AR"=rnorm(n=1,mean=0,sd=.1),
               "beta3_AR"=rnorm(n=1,mean=0,sd=1),
               "beta1_EN"=rnorm(n=1,mean=0,sd=1),
               "beta2_EN"=rnorm(n=1,mean=0,sd=.1),
               "beta3_EN"=rnorm(n=1,mean=0,sd=1),
               "prec_AR"=runif(n=1,0,10),
               "prec_EN"=runif(n=1,0,10),
               "rho"=1,
               "kappa"=runif(n=1,0,1)))
} 

samples <- runMCMC(ChospitalMCMC,  niter =150000,nchains = 3,nburnin=50000
                   ,samplesAsCodaMCMC = TRUE,thin=10,inits = initsFunction)


##test convergence
#should all be less than 1.05
gelman.diag(samples[,c("beta0_AR","beta1_AR",
                       "beta2_AR","beta3_AR",
                       "r","prec_AR",
                       "beta0_EN","prec_EN","rho","kappa",
                       "beta1_EN",
                       "beta2_EN","beta3_EN")])

#should be greater than 1000
min(effectiveSize(samples[,c("beta0_AR","beta1_AR",
                       "beta2_AR","beta3_AR",
                       "r","prec_AR",
                       "beta0_EN","prec_EN","rho","kappa",
                       "beta1_EN",
                       "beta2_EN","beta3_EN")]))


#can examine some traceplots
plot(samples[,c("beta0_AR","beta1_AR","rho")])
plot(samples[,c("beta2_AR","beta3_AR")])
plot(samples[,c("beta2_EN","beta3_EN")])
plot(samples[,c("beta1_EN","beta1_AR")])
plot(samples[,c("r","prec_AR")])
plot(samples[,c("beta0_EN","prec_EN")])
plot(samples[,c("beta0_AR","prec_AR")])
plot(samples[,c("kappa","kappa_prior")])
plot(samples[,c("weight[1]","weight[2]")])
plot(samples[,c("weight_prior[1]","weight_prior[2]")])


#now can calculate the WAIC from Table 1
mmobility_matrix <- mean(mobility_matrix)
N_matrix <- matrix(rep(0,30*30),nrow=30,ncol=30)
for(i in 1:30){
  N_matrix[i,near_nei_matrix[i,]] <- 1
}
diag(N_matrix) <- 1
samps <- data.frame(rbind(samples[[1]],samples[[2]],samples[[3]]))
lppd <- 0 
pwaic <- 0

for(i in 1:N){
  print(i)
  for(t in 3:mT){
    
    #calculate nei weighted sum
    y_nei_sumitm1 <- rep(0,30000)
    neii <- which(N_matrix[i,]==1)
    for(j in neii){
      y_nei_sumitm1 <- y_nei_sumitm1 + psi[j,t-1]*samps[,paste0("W_norm.",i,"..",j,".")]
      
    }
    y_nei_sumitm2 <- rep(0,30000)
    neii <- which(N_matrix[i,]==1)
    for(j in neii){
      y_nei_sumitm2 <- y_nei_sumitm2 + psi[j,t-2]*samps[,paste0("W_norm.",i,"..",j,".")]
      
    }
    
    
    lambda <- exp(samps[,paste0("b_AR.",i,".")]+
                    samps[,"beta2_AR"]*(mobility_matrix[i,t-1]-mmobility_matrix)+
                    samps[,"beta3_AR"]*new_variant[t])*(y_nei_sumitm1*samps[,"weight.1."]+
                                                          y_nei_sumitm2*samps[,"weight.2."])+
      exp(samps[,paste0("b_EN.",i,".")]+
            samps[,"beta2_EN"]*(mobility_matrix[i,t-1]-mmobility_matrix)+
            samps[,"beta3_EN"]*new_variant[t])
    
    
    
    r<- samps[,"r"]
    
    
    lppd <- lppd + log(mean(dnbinom(x = admissions[i,t],prob = r/(r+lambda),size = r)))
    pwaic <- pwaic + var(log(dnbinom(x = admissions[i,t],prob = r/(r+lambda),size = r)))
  }
}

waic <- -2*(lppd-pwaic)
#should be 17,845

#fitted values
fitted <- array(dim = c(N,mT,30000))

for(i in 1:N){
  print(i)
  for(t in 3:mT){
    
    #calculate nei weighted sum
    y_nei_sumitm1 <- rep(0,30000)
    neii <- which(N_matrix[i,]==1)
    for(j in neii){
      y_nei_sumitm1 <- y_nei_sumitm1 + psi[j,t-1]*samps[,paste0("W_norm.",i,"..",j,".")]
      
    }
    y_nei_sumitm2 <- rep(0,30000)
    neii <- which(N_matrix[i,]==1)
    for(j in neii){
      y_nei_sumitm2 <- y_nei_sumitm2 + psi[j,t-2]*samps[,paste0("W_norm.",i,"..",j,".")]
      
    }
    
    
    lambda <- exp(samps[,paste0("b_AR.",i,".")]+
                    samps[,"beta2_AR"]*(mobility_matrix[i,t-1]-mmobility_matrix)+
                    samps[,"beta3_AR"]*new_variant[t])*(y_nei_sumitm1*samps[,"weight.1."]+
                                                          y_nei_sumitm2*samps[,"weight.2."])+
      exp(samps[,paste0("b_EN.",i,".")]+
            samps[,"beta2_EN"]*(mobility_matrix[i,t-1]-mmobility_matrix)+
            samps[,"beta3_EN"]*new_variant[t])
    
    
    
    r<- samps[,"r"]
    
    
    
    fitted[i,t,] <- rnbinom(n = 30000,
                            prob=r/(r+lambda),size=r)
    
  }
}

#can plot the fit
for(i in 1:30){
  
  med = apply(fitted[i,,],MARGIN = 1,mean)
  upper= apply(fitted[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.975),na.rm=TRUE))
  lower= apply(fitted[i,,],MARGIN = 1,function(x) quantile(x,probs=c(.025),na.rm=TRUE))
  
  
  plot(1:mT,admissions[i,],type = "p",main=paste0(i),ylim=c(0,max(c(upper),na.rm = TRUE)+.25*mean(c(med),na.rm = TRUE)),
       cex=.5)
  
  
  lines(1:mT,med,col="red")
  lines(1:mT,lower,col="red",lty=2)
  lines(1:mT,upper,col="red",lty=2)
  
}



