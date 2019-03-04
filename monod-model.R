

library(tidyverse)
library(cowplot)
library(simecol)


data <- read_csv("data-raw/population4-nitrate.csv") %>% 
	select(days, RFU, nitrate_concentration, population, well_plate) %>% 
	rename(starting_nitrate = nitrate_concentration) %>% 
	mutate(nitrate = starting_nitrate - RFU/20) %>% 
	rename(N = RFU,
		   R = nitrate) %>% 
	filter(starting_nitrate == 1000, well_plate == "B03_29") 

data %>%
	ggplot(aes(x = days, y = N)) + geom_point() +
	geom_point(data = data, aes(x = days, y = R), color = "blue") 

Parameters <- c(r = 3, k = 15, m = 0.1)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.01, K = 0.1, m = 0)
UpperBound <- c(r = 3, K = 500, m = 1) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound


monod_model <- new("odeModel",
			   main = function (time, init, parms) {
			   	with(as.list(c(init, parms)), {
			   		dN <-  r * (R / (R + k)) - m
			   		dR <-  (-N + m * N) / 1
			   		list(c(dN, dR))
			   	})
			   },
			   parms = Parameters,
			   times = c(from = 0, to = 3, by = 0.1), # the time interval over which the model will be simulated.
			   init = c(N = 15, R = 999.25),
			   solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
)


# These vectors simply contain strings, which are used to facilitate parameter
# assignment in the below code. While seemingly clumsy, it appears to be a
# necessary step.
fittedparms <- c("r", "k", "m") # for assigning fitted parameter values to fittedCRmodel

## Model Fitting Function ##

# The following function is intended to be used with map_df() on the nested
# dataframe called "controldata". It takes a single dataframe of various
# observations for a control replicate, and outputs a dataframe consisting of
# the replicate ID, the Phosphorus treatment, the temperature, and the
# parameter estimates for r and K. It can also be used to output parameter
# values for a single replicate. To do this call controlfit(controldata[['X']],
# where "X" is the replicate's ID number.

monod_fit <- function(data){
	
	init(monod_model) <- c(N = data$N[1], R = data$R[1]) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, N) # The Y values of the observed data points we are fitting our model to
	parms(monod_model) <- Parameters
	# parms(monod_model)[TempName] <- data$temp[1] # Set the temperature parameter in CRmodel to whatever our control replicate used.
	
	# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
	# squared differences between the experimental data and our modelled data. This
	# is fairly standard, although alternatives do exist.
	
	# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
	# "lower" is a vector containing the lower bound constraints
	# for the parameter values. This may need tweaking.
	
	fitted_monod_model <- fitOdeModel(monod_model, whichpar = fittedparms, obstime, yobs,
								 debuglevel = 0, fn = ssqOdeModel,
								 method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
								 control = list(trace = T)
	)
	
	# Here we create vectors to be used to output a dataframe of
	# the replicates' ID, Phosphorus level, temperature, and the
	# fitted parameters. "truer" and "trueK" are the fitted
	# parameters, but scaled using the appropriate arrhenius
	# transform.	
	population <- data$population[1]
	r <- coef(fitted_monod_model)[1]
	k <- coef(fitted_monod_model)[2]
	m <- coef(fitted_monod_model)[3]
	ID <- data$well_plate[1]
	# population <- tail(parms(CRmodel), n=1)
	output <- data.frame(ID, population, r, k, m)
	return(output)
}

# Here we fit the values for r and K, using the 12C replicates, and then
# prepare these data to be used to fit the activation energies further below in
# the code.

# Fit r and K for 12C replicates
fits <- map_df(data, monod_fit)

summary(fitted_monod_model)

rKfulldata <- map_df(TwelveFullPdata, rKfit)


# These vectors simply contain strings, which are used to facilitate parameter
# assignment in the below code. While seemingly clumsy, it appears to be a
# necessary step.
TempName <- c("temp") # for assigning temperature values to CRmodel
fittedparms <- c("r", "K") # for assigning fitted parameter values to fittedCRmodel

## Model Fitting Function ##

# The following function is intended to be used with map_df() on the nested
# dataframe called "controldata". It takes a single dataframe of various
# observations for a control replicate, and outputs a dataframe consisting of
# the replicate ID, the Phosphorus treatment, the temperature, and the
# parameter estimates for r and K. It can also be used to output parameter
# values for a single replicate. To do this call controlfit(controldata[['X']],
# where "X" is the replicate's ID number.

rKfit <- function(data){
	
	init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
	parms(CRmodel) <- Parameters
	parms(CRmodel)[TempName] <- data$temp[1] # Set the temperature parameter in CRmodel to whatever our control replicate used.
	
	# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
	# squared differences between the experimental data and our modelled data. This
	# is fairly standard, although alternatives do exist.
	
	# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
	# "lower" is a vector containing the lower bound constraints
	# for the parameter values. This may need tweaking.
	
	fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
								 debuglevel = 0, fn = ssqOdeModel,
								 method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
								 control = list(trace = T)
	)
	
	# Here we create vectors to be used to output a dataframe of
	# the replicates' ID, Phosphorus level, temperature, and the
	# fitted parameters. "truer" and "trueK" are the fitted
	# parameters, but scaled using the appropriate arrhenius
	# transform.	
	population <- data$population[1]
	r <- coef(fittedCRmodel)[1]
	K <- coef(fittedCRmodel)[2]
	ID <- data$ID[1]
	starting_nitrate <- tail(parms(CRmodel), n=1)
	output <- data.frame(ID, population, starting_nitrate, r, K)
	return(output)
}

# Here we fit the values for r and K, using the 12C replicates, and then
# prepare these data to be used to fit the activation energies further below in
# the code.

# Fit r and K for 12C replicates
rKdefdata <- map_df(TwelveDefPdata, rKfit)
rKfulldata <- map_df(TwelveFullPdata, rKfit)



pfit <- function(data){
	
	day_three_cell_volume <- data$volume_cell[1]
	day_zero_cell_concentration <- 10 ^ 5
	initial_algal_biovolume <- day_three_cell_volume * day_zero_cell_concentration
	data <- add_row(data, H = 10, P = initial_algal_biovolume, days = 0, .before = 1)		
	
	temp <- data$temperature[2]
	model <- CRmodel
	init(model) <- c(P = data$P[1], H = 10) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, P, H) # The Y values of the observed data points we are fitting our model to
	
	# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
	# squared differences between the experimental data and our modelled data. This
	# is fairly standard, although alternatives do exist.
	
	# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
	# "lower" is a vector containing the lower bound constraints
	# for the parameter values. This may need tweaking.
	
	fittedmodel <- fitOdeModel(model, whichpar = FittedParameters, obstime, yobs,
							   debuglevel = 0, fn = ssqOdeModel,
							   method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
							   control = list(trace = TRUE),
							   rtol = 1e-9,
							   atol = 1e-9
	)
	

	
	# Here we create vectors to be used to output a dataframe of
	# the replicates' ID, Phosphorus level, temperature, and the
	# fitted parameters. 
	
	ID <- data$unique_ID[2]
	population <- data$population[2]
	r <- coef(fittedmodel)["r"]
	K <- coef(fittedmodel)["K"]
	a <- coef(fittedmodel)["a"]
	b <- coef(fittedmodel)["b"]
	eps <- coef(fittedmodel)["eps"]
	m <- coef(fittedmodel)["m"]
	
	simmodel <- model
	# set model parameters to fitted values and simulate again
	parms(simmodel)[FittedParameters] <- coef(fittedmodel)
	simdata <- out(sim(simmodel, rtol = 1e-9, atol = 1e-9))
	
	finalsimulatedP <- last(simdata$P)
	finalsimulatedH <- last(simdata$H)
	finalobservedP <- last(yobs$P)
	finalobservedH <- last(yobs$H)
	
	output <- data.frame(ID, Phosphorus, temp, r, K, a, b, eps, m,
						 finalsimulatedP,
						 finalsimulatedH,
						 finalobservedP,
						 finalobservedH)
	
	return(output)
}


