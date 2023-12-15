#will create the simulated data set used in Section 4 of the main text
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

#initializes the simulated counts y and states S
ysim <- matrix(nrow=30,ncol=120)
Ssim <- matrix(nrow=30,ncol=120)

#need to make clusters
#5 clusters of 6 areas, so each area has 5 neighbors
cluster <- rep(NA,30)
cluster[1:6] <- 1
cluster[7:12] <- 2
cluster[13:18] <- 3
cluster[19:24] <- 4
cluster[25:30] <- 5

##generates the states
#2=endemic and 3=outbreak
#For every area, the start of each outbreak was
#randomized to occur within the first 4 weeks of the corresponding outbreak period.
#this was done to create spatial synchronization in the outbreaks within a cluster
for(i in 1:30){
  
  Ssim[i,1:15] <- 2
  start_one <- sample(1:4,size=1)
  if(start_one==1){
    Ssim[i,16:30] <- 3
  }else{
    Ssim[i,16:(16+(start_one-2))] <- 2
    Ssim[i,(16+(start_one-1)):30] <- 3
  }
  
  Ssim[,31:45] <- 2
  start_two <- sample(1:4,size=1)
  if(start_two==1){
    Ssim[i,46:60] <- 3
  }else{
    Ssim[i,46:(46+(start_two-2))] <- 2
    Ssim[i,(46+(start_two-1)):60] <- 3
  }
  
  Ssim[,61:75] <- 2
  start_three <- sample(1:4,size=1)
  if(start_three==1){
    Ssim[i,76:90] <- 3
  }else{
    Ssim[i,76:(76+(start_three-2))] <- 2
    Ssim[i,(76+(start_three-1)):90] <- 3
  }
  
  Ssim[,91:105] <- 2
  start_four <- sample(1:4,size=1)
  if(start_four==1){
    Ssim[i,106:120] <- 3
  }else{
    Ssim[i,106:(106+(start_four-2))] <- 2
    Ssim[i,(106+(start_four-1)):120] <- 3
  }
  
  
}

#now enter the absence periods
#1=absence
#we assumed a 40 percent chance each endemic period had
#an absence period inserted into the middle
for(i in 1:30){
  for(k in c(1,3,5,7)){
    if(runif(n=1)<.4){
      Ssim[i,((k-1)*15+5):(k*15-4)] <- 1
    }
  }
}

#make sure no NAs
sum(is.na(Ssim))

#beds covariate
standard_cap <- c(772,282,359,380,321,253,531,320,240,376,487,288,421,244,419,104,
                  257,261,564,256,191,544,274,371,373,341,277,293,324,467)/100

#now we generate y conditional on the states
#we assume previously there are were no cases
ysim[,1] <- rnbinom(n=30,size=10,
                    mu=exp(0+(standard_cap-mean(standard_cap))*.1))
for(i in 1:30){
  for(t in 2:120){
    if(Ssim[i,t]==1){
      ysim[i,t] <- 0
    }else if(Ssim[i,t]==2) {
      ysim[i,t] <- rnbinom(n=1,size=10,
                           mu=exp(0+(standard_cap[i]-mean(standard_cap))*.1)*((ysim[i,t-1]+1)^(.5))
      )
    }else if(Ssim[i,t]==3) {
      ysim[i,t] <- rnbinom(n=1,size=10,
                           mu=exp(.75+(standard_cap[i]-mean(standard_cap))*.05)*((ysim[i,t-1]+1)^(.75))
      )
    }
  }
}

#a table of the states
#should make sure there are not too few of one state
#or model fitting may be unstable
table(Ssim)

#can plot the simulated data in some areas
plot(ysim[1,])
lines(ysim[1,])
plot(Ssim[1,])
plot(ysim[6,])
lines(ysim[6,])
plot(Ssim[6,])



