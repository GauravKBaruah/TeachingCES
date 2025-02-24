rm(list=ls())

#R script for modelling simple eco-evolutionary dynamics of a population

#Ecological model - no evolution of the phenotype, but the phenotype of the insect impacts growth rate
# Population dynamics only
b0<-1.5
d0<-0.5
w<-2  #width of the fitness curve
T<-20 #optimum Temperature
K<-200 #carrying capacity
u<- 20 # mean average trait value (the average optimum temperature at which fitness is the highest)
s<-0.1 #genetic variance. Here we assume environmental variance to be zero.
tmax<-500 #100 time points
dt<-0.05
N<-numeric(tmax)
N[1]<-100

for(t in 1:tmax){
  b<- (b0*w/sqrt(w^2+2*s))*exp(-(u-T)^2/(2*s+w^2))  #average rate of growth
  N[t+1] <- N[t] + N[t]*(b- d0- N[t]/K)*dt #euler approximation. This is being done in a simpler way. Euler approximation
  #is not accurate and one should switch to runge-kutta 4 or solvers for better accuracy.
  
  if(N[t+1]<0){N[t+1]<-0} # making sure abundance does drop below 0
}
plot(0:tmax,N, ylab="Abundance",xlab="Time",typ='l',lwd=3)




### What happens when Temperature is now 22 C ?


#Ecological model - no evolution of the phenotype, but the phenotype of the insect impacts growth rate
# Population dynamics only
b0<-1.5
d0<-0.5
w<-2 #width of the fitness curve
T<-21 #optimum Temperature
K<-200 #carrying capacity
u<- 20 # mean average trait value (the average optimum temperature at which fitness is the highest)
s<-0.1 #genetic variance. Here we assume environmental variance to be zero.
tmax<-500 #100 time points
dt<-0.05
N<-numeric(tmax)
N[1]<-190

for(t in 1:tmax){
  b<- (b0*w/sqrt(w^2+2*s))*exp(-(u-T)^2/(2*s+w^2))
  N[t+1] <- N[t] + N[t]*(b- d0-N[t]/K)*dt #euler approximation. This is being done in a simpler way. Euler approximation
  #is not accurate and one should switch to runge-kutta 4 or solvers for better accuracy.
  if(N[t+1]<0){N[t+1]<-0}
  
}
lines(0:tmax,N, ylab="Abundance",xlab="Time",typ='l', lwd=3,col="firebrick")
legend("topright", lty=c(1,1),col=c("black","firebrick"), c("no mismatch", "mismatch"))



############# Eco-Evolutionary dynamics ##############

# Here not only the population changes over time but also the mean trait evolves in response to changes

b0<-1.5
w<-2 #width of the fitness curve
T<-20 #optimum temperature
K<-200 #carrying capacity
s<-0.1 #genetic variance. Here we assume environmental variance to be zero. So Vp=Vg
tmax<-500 #100 time points
dt<-0.1
N<-numeric(tmax)
u<-numeric(tmax)
u[1]<-20# mean average trait value (the average optimum temperature at which fitness is the highest)
N[1]<-100

for(t in 1:tmax){
  b<- (b0*w/sqrt(w^2+2*s))*exp(-(u[t]-T)^2/(2*s+w^2))  #average rate of growth #average rate of growth
  N[t+1] <- N[t] + N[t]*(b- d0-N[t]/K)*dt #euler approximation. This is being done in a simpler way. Euler approximation
  #is not accurate and one should switch to runge-kutta 4 or solvers for better accuracy.
  dbdu<- -2*w*(u[t]-T)*exp(-(u[t]- T)^2/(2*s+w^2))/(2*s+w^2)^3/2
  u[t+1]<- u[t]+ (s*dbdu)*dt   #V_g*dr/du
  
}
par(mfrow=c(1,2))
plot(0:tmax,N, ylab="Abundance",xlab="Time",typ='l', lwd=3,col="firebrick")
plot(0:tmax,u, ylab="Mean trait",xlab="Time",typ='l', lwd=3,col="firebrick", ylim = c(19,21))




### 0.5 degree mismatch
b0<-1.5
w<-2 #width of the fitness curve
T<-22 #average local Temperature
K<-500 #carrying capacity
s<-1 #genetic variance. Here we assume environmental variance to be zero.
tmax<-500 #100 time points
dt<-0.1
N<-numeric(tmax)
u<-numeric(tmax)
u[1]<-20# mean average trait value (the average optimum temperature at which fitness is the highest)
N[1]<-193

for(t in 1:tmax){
  b<- (b0*w/sqrt(w^2+2*s))*exp(-(u[t]-T)^2/(2*s+w^2))  #average rate of growth #average rate of growth
  N[t+1] <- N[t] + N[t]*(b- d0-N[t]/K)*dt #euler approximation. This is being done in a simpler way. Euler approximation
  #is not accurate and one should switch to runge-kutta 4 or solvers for better accuracy.
  dbdu<- -2*w*(u[t]-T)*exp(-(u[t]- T)^2/(2*s+w^2))/(2*s+w^2)^3/2
  u[t+1]<- u[t]+ (s*dbdu)*dt   #V_g*dr/du
  
}
par(mfrow=c(1,2))
plot(0:tmax,N, ylab="Abundance",xlab="Time",typ='l', lwd=3,col="firebrick")
plot(0:tmax,u, ylab="Mean trait",xlab="Time",typ='l', lwd=3,col="firebrick")




############ Some random stochasticity in dynamics of population ###################




### 0.5 degree mismatch
b0<-1.5
w<-2 #width of the fitness curve
T<-21.5 #average local Temperature
K<-500 #carrying capacity
s<-1 #genetic variance. Here we assume environmental variance to be zero.
tmax<-500 #100 time points
dt<-0.1
N<-numeric(tmax)
u<-numeric(tmax)
u[1]<-20# mean average trait value (the average optimum temperature at which fitness is the highest)
N[1]<-100

noise_env<-rnorm(tmax,0,0.25)
for(t in 1:tmax){
  b<- (b0*w/sqrt(w^2+2*s))*exp(-(u[t]-T)^2/(2*s+w^2))  #average rate of growth #average rate of growth
  N[t+1] <- N[t] + N[t]*(b- d0-N[t]/K)*dt + N[t]*noise_env[t]*dt#euler approximation. This is being done in a simpler way. Euler approximation
  #is not accurate and one should switch to runge-kutta 4 or solvers for better accuracy.
  dbdu<- -2*w*(u[t]-T)*exp(-(u[t]- T)^2/(2*s+w^2))/(2*s+w^2)^3/2
  u[t+1]<- u[t]+ (s*dbdu)*dt   #V_g*dr/du
  
}
par(mfrow=c(1,2))
plot(0:tmax,N, ylab="Abundance",xlab="Time",typ='l', lwd=3,col="firebrick")
plot(0:tmax,u, ylab="Mean trait",xlab="Time",typ='l', lwd=3,col="firebrick")



######################################################################### HIGH MISMATCH HIGH GENETIC VARIANCE ###########


### 2.5 degree mismatch and high genetic variance
b0<-1.5
w<-2 #width of the fitness curve
T<-22.5 #average local Temperature
K<-500 #carrying capacity
s<-1#genetic variance. Here we assume environmental variance to be zero and genetic variance to be high
tmax<-1000 #100 time points
dt<-0.1
N<-numeric(tmax)
u<-numeric(tmax)
u[1]<-20# mean average trait value (the average optimum temperature at which fitness is the highest)
N[1]<-500

noise_env<-rnorm(tmax,0,0.25)
for(t in 1:tmax){
  b<- (b0*w/sqrt(w^2+2*s))*exp(-(u[t]-T)^2/(2*s+w^2))  #average rate of growth #average rate of growth
  N[t+1] <- N[t] + N[t]*(b- d0-N[t]/K)*dt + N[t]*noise_env[t]*dt#euler approximation. This is being done in a simpler way. Euler approximation
  #is not accurate and one should switch to runge-kutta 4 or solvers for better accuracy.
  dbdu<- -2*w*(u[t]-T)*exp(-(u[t]- T)^2/(2*s+w^2))/(2*s+w^2)^3/2
  u[t+1]<- u[t]+ (s*dbdu)*dt  #V_g*dr/du
  
}
par(mfrow=c(1,2))
plot(0:tmax,N, ylab="Abundance",xlab="Time",typ='l', lwd=3,col="firebrick")
plot(0:tmax,u, ylab="Mean trait",xlab="Time",typ='l', lwd=3,col="firebrick")




######################################################################### LOW GENETIC VARIANCE HIGH MISMATCh###############


### 2.5 degree mismatch and low genetic variance
b0<-1.5
w<-2 #width of the fitness curve
T<-22.5 #average local Temperature
K<-500 #carrying capacity
s<-0.2 #low genetic variance. Here we assume environmental variance to be zero and genetic variance to be high
tmax<-1000 #100 time points
dt<-0.1
N<-numeric(tmax)
u<-numeric(tmax)
u[1]<-20# mean average trait value (the average optimum temperature at which fitness is the highest)
N[1]<-500

noise_env<-rnorm(tmax,0,0.25)
for(t in 1:tmax){
  b<- (b0*w/sqrt(w^2+2*s))*exp(-(u[t]-T)^2/(2*s+w^2))  #average rate of growth #average rate of growth
  N[t+1] <- N[t] + N[t]*(b- d0-N[t]/K)*dt + N[t]*noise_env[t]*dt#euler approximation. This is being done in a simpler way. Euler approximation
  #is not accurate and one should switch to runge-kutta 4 or solvers for better accuracy.
  dbdu<- -2*w*(u[t]-T)*exp(-(u[t]- T)^2/(2*s+w^2))/(2*s+w^2)^3/2
  u[t+1]<- u[t]+ (s*dbdu)*dt  #V_g*dr/du
  
}
par(mfrow=c(1,2))
plot(0:tmax,N, ylab="Abundance",xlab="Time",typ='l', lwd=3,col="firebrick")
plot(0:tmax,u, ylab="Mean trait",xlab="Time",typ='l', lwd=3,col="firebrick")

