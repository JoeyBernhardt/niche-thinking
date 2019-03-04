

library(tidyverse)
library(cowplot)
library(simecol)


data <- read_csv("data-raw/population4-nitrate.csv") %>% 
	select(days, RFU, nitrate_concentration, population, well_plate) %>% 
	rename(starting_nitrate = nitrate_concentration) %>% 
	mutate(nitrate = starting_nitrate - RFU/20) %>% 
	rename(N = RFU,
		   R = nitrate) %>% 
	filter(starting_nitrate == 1000, well_plate == "F03_29") 

data <- read_csv("data-processed/pop4_b03.csv") ## manually added a bit of error here

data %>%
	ggplot(aes(x = days, y = N)) + geom_point() +
	geom_point(data = data, aes(x = days, y = R), color = "blue") 



# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.01, k = 0.1)
UpperBound <- c(r = 3, k = 500) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound




chemostat


Parameters <- c(r = 1.5, k = 0.1)
monod_model <- new("odeModel",
				   main = function (time, init, parms) {
				   	with(as.list(c(init, parms)), {
				   		dN <-  r * (R / (R + k)) - 0.1
				   		dR <-  (-N + 0.1 * N) / 1
				   		list(c(dN, dR))
				   	})
				   },
				   parms = Parameters,
				   times = c(from = 0, to = 3, by = 0.1), # the time interval over which the model will be simulated.
				   init = c(N = 15, R = 999.25),
				   solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
)

lv <- sim(monod_model)
plot(lv)

# These vectors simply contain strings, which are used to facilitate parameter
# assignment in the below code. While seemingly clumsy, it appears to be a
# necessary step.
fittedparms <- c("r", "k") # for assigning fitted parameter values to fittedCRmodel

monod_fit <- function(data){
	
	init(monod_model) <- c(N = data$N[1], R = data$R[1]) # Set initial model conditions 
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, N, R) # The Y values of the observed data points we are fitting our model to
	parms(monod_model) <- Parameters

	fitted_monod_model <- fitOdeModel(monod_model, whichpar = fittedparms, obstime, yobs,
									  debuglevel = 0, fn = ssqOdeModel,
									  method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
									  control = list(trace = TRUE),
									  rtol = 1e-9,
									  atol = 1e-9
	)
	
	population <- data$population[1]
	r <- coef(fitted_monod_model)[1]
	k <- coef(fitted_monod_model)[2]
	# m <- coef(fitted_monod_model)[3]
	ID <- data$well_plate[1]
	# population <- tail(parms(CRmodel), n=1)
	output <- data.frame(ID, population, r, k)
	return(output)
}

simmodel <- monod_model
parms(simmodel)[FittedParameters] <- coef(fitted_monod_model)
simdata <- out(sim(simmodel, rtol = 1e-9, atol = 1e-9))

parms(monod_model)
# Here we fit the values for r and K, using the 12C replicates, and then
# prepare these data to be used to fit the activation energies further below in
# the code.

# Fit r and K for 12C replicates
fits <- map_df(data, monod_fit)
