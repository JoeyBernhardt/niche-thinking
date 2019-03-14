


egdat<-data.frame(time=c(17.9284722,9.9145833,13.9541667,9.0125000,12.9284722,16.9291667,0.9722222,11.0013889,0.0000000,8.0131944,11.9333333,18.9062500,4.0229167,3.0333333,7.0062500,5.9847222,20.9340278,19.9395833,4.8868056,14.9965278,16.0180556,1.8909722),
				  log.fluor= exp(c(5.112714,4.948780,5.096441,4.588206,5.092295,5.119763,1.449947,5.110046,0.781819,4.255234,4.998089,5.055782,2.663951,2.328480,3.899183,3.413726,4.966042,4.953977,2.900249,5.226657,5.121840,1.808143)))

par(mfrow=c(1,1))
plot(log.fluor~time,data=egdat)

### r - alpha formulation and fitting


parameters<-c(r = 1, alpha = 0.001)
state<-c(N=exp(0.781819))
logistic <- function(t,state,parameters){
	with(as.list(c(state,parameters)),{
		dN <- (r - alpha*N)*N
		
		# return the rate of change
		list(c(dN))
	})	# end of with(as.list...
}

times<-seq(0,20,0.001)

out<-ode(y=state,times=times,func= logistic,parms=parameters)
out<-data.frame(out)

# plot solution
par(mfrow=c(2,1))
plot(N~time,data=out,type='l',col='red')
plot(log10(N)~time,data=out,type='l',col='red')


parameters<-c(r = 0.86, alpha = 0.005)
state<-c(N=100)
# Wrap Monod model ODE in a function that's easier to use with bbmle/mle2
r_alpha_ode<-function(time.vals, N0, r, alpha,verbose=F,maxtime){
	
	parameters<-c(r=r, alpha = alpha)
	state<-c(N=N0)
	
	if(verbose){
		print(c(parameters,state))
	}
	
	# differential equation set up

	r_alpha <- function(t,state,parameters){
		with(as.list(c(state,parameters)),{
			dN <- (r - alpha*N)*N
			
			# return the rate of change
			list(c(dN))
		})	# end of with(as.list...
	}
	times<-unique(sort(c(seq(0,maxtime,0.005),time.vals)))
	
	out<-ode(y=state,times=times,func = r_alpha,parms=parameters)
	out<-data.frame(out)
	
	# take values from ode solution as predicted values for the given observation times
	time.df<-data.frame(id=seq(1,length(time.vals)),time=time.vals)
	out2<-merge(time.df,out,by='time',all.x=T)
	out2<-out2[order(out2$id),]
	
	preds<-out2$N
	preds<-log(preds)
	
	return(preds)
}


### Fit whole data set: ###

parameters<-c(r = 1, alpha = 0.001)
state<-c(N=exp(0.781819))

maxtime<-max(egdat$time)

N0.guess<- min(egdat$log.fluor)
exp(0.781819)


### Run fit (slow)

# With my first set of starting guesses, the fit doesn't quite converge
fit_r_alpha_1<-mle2(log.fluor ~ dnorm(mean = r_alpha_ode(time, N0, r, alpha, maxtime=maxtime,verbose = T), sd=exp(s)),
					start=list(N0= N0.guess, r=1.5, alpha = 0.001, s =-0.1), data=egdat)


# Note that parameters N0, m, R0, and K are all constrained to be positive by wrapping their values in the exp() function; as a result, the optimization algorithm can't force their values to be negative.

summary(fit_r_alpha_1)
# grab the final values
# r           m          R0           K           N 
# 0.6052132   0.1359362 164.6618541  14.3820060   2.5874721 

# run again with updated starting guesses (works)
fit.monod.2<-mle2(log.fluor~dnorm(mean=monod.ode(time,exp(N0),r,exp(m),exp(R0),exp(K),maxtime=maxtime),sd=exp(s)),
				  start=list(N0=log(2.5874721),r=0.6052132,m=log(0.1359362),R0=log(164.6618541),K=log(14.3820060),s=-0.1),data=egdat)

summary(fit.monod.2)

# Coefficients:
#      Estimate Std. Error  z value  Pr(z)    
#   N0  0.95058    0.04998  19.0193 <2e-16 ***
#   r   0.61695    1.49972   0.4114 0.6808    
#   m  -1.91446   10.18222  -0.1880 0.8509    
#   R0  5.10564    0.21080  24.2202 <2e-16 ***
#   K   2.65149    2.65907   0.9972 0.3187    
#   s  -2.56094    0.15074 -16.9891 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: -50.24267

# extract and save coefficients
cfs<-unlist(coef(fit.monod.2))
attr(cfs,"names")<-NULL

N0.est<-exp(cfs[1])
r.est<-cfs[2]
m.est<-exp(cfs[3])
R0.est<-exp(cfs[4])
K.est<-exp(cfs[5])


### Plot the results:

# Set up parameters and initial conditions
parameters<-c(r=r.est,K=K.est,Y=1,m=m.est,R0=R0.est)
state<-c(N=N0.est)

Monod<-function(t,state,parameters){
	with(as.list(c(state,parameters)),{
		dN <- (r*((R0-Y*N)/((R0-Y*N)+K))-m)*N	
		
		# return the rate of change
		list(c(dN))
	})	# end of with(as.list...
}
times<-seq(0,22,0.001)
out<-ode(y=state,times=times,func=Monod,parms=parameters)
out<-data.frame(out)

# Plot results:
plot(log.fluor~time,data=egdat)
lines(log(N)~time,data=out,col='blue')

