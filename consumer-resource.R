### growth functions
library(deSolve)
library(tidyverse)
library(cowplot)




R1_input = 10
R2_input = 12
D1_input = 1
D2_input = 1
N1_input = 40
N2_input = 40
f11_input = 4
f21_input = 0.4
f12_input = 3
f22_input = 0.4
a12_input = 0.2
a22_input = 0.3
a11_input = 0.2
a21_input = 0.4
S_input = 10
c1_input = 0.5
c2_input = 0.5

Rstar1 <- D1_input/(f11_input*a11_input)
ZNGI_1_slope <- (-f11_input*a11_input)/(f12_input*a12_input)
ZNGI_1_intercept <- (D1_input - f11_input*a11_input*Rstar1)/(f12_input*a12_input)


Rstar2 <- D2_input/(f21_input*a21_input)
ZNGI_2_slope <- (-f21_input*a21_input)/(f22_input*a22_input)
ZNGI_2_intercept <- (D2_input - f21_input*a21_input*Rstar2)/(f22_input*a22_input)

ggplot() +
	geom_abline(slope = ZNGI_1_slope, intercept = ZNGI_1_intercept, color = "green") +
	geom_abline(slope = ZNGI_2_slope, intercept = ZNGI_2_intercept, color = "blue") +
	ylim(0, 10) +xlim(0, 10) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + coord_cartesian()



# logistic <- function(R1_input = 10, R2_input = 20, D1_input = 0,  D2_input = 0,
# 					 N1_input = 40, N2_input = 40, f11_input = 0.4, f21_input = 0.3,
# 					 f12_input = 0.3, f22_input = 0.3, 
# 					 a12_input = 0.7, a22_input = 0.1,
# 					 a11_input = 0.2, a21_input = 0.5,
# 					 S_input = 5, c1_input = 0.2, c2_input = 0.2){
	
	dNRdt <-  function(t, state, parameters) {
		with(
			as.list(c(state,parameters)),{
				dN1 = N1 * (f11*a11*R1 + f12*a12*R2 - D1) 
				dN2 = N2 * (f21*a21*R1 + f22*a22*R2 - D2)  #algae
				dR1 = c1 * (S - R1) - (f11*N1*R1 + f21*N2*R1) 
				dR2 = c2 * (S - R2) - (f12*N1*R2 + f22*N2*R2) #resource
				return(list(c(dN1, dN2, dR1, dR2)))			
			}
		)
	}
	
	
	parameters <-  c(R1 = R1_input, R2= R2_input, D1 = D1_input, D2 = D2_input,
					 N1= N1_input, N2 = N2_input, f11 = f11_input, f21 = f21_input,
					 f12 = f12_input, f22 = f22_input, 
					 a12 = a12_input, a22 = a22_input,
					 a11 = a11_input, a21 = a21_input,
					 S = S_input, c1 = c1_input, c2 = c2_input)
	state <-  c(N1 = N1_input, N2 = N2_input, R1 = R1_input,  R2 = R2_input) ## initial conditions
	time <-  seq(from = 0, to = 50, by = 0.01) 
	
	## solve numerically
	out <- ode(y = state, times = time, func = dNRdt, parms = parameters)
	
	## tidy up the output and plot
	out.df <- as.data.frame(out) %>%
		rename(algae1 = N1,
			   algae2 = N2,
			   resource1 = R1,
			   resource2 = R2) %>% 
		gather(key = type, value = value, algae1, algae2, resource1, resource2)
	# Return to the main program, passing back the timeseries of P   
	# return(out.df)
# }




out.df %>% 
	# filter(type == "resource2") %>% 
	ggplot(aes( x= time, y = value, color = type)) + geom_line()


