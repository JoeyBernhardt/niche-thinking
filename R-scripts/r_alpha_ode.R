library(simecol)
library(tidyverse)
library(cowplot)


### notes march 20 2019 -- keep getting the alpha estimates hitting the lower bounds ugh. not sure what to do next. 


Parameters <- c(r = 1, alpha = 0.01) ## initial parameter guesses
CRmodel <- new("odeModel",
			   main = function (time, init, parms) {
			   	with(as.list(c(init, parms)), {
			   		dN <-  (r - alpha*N)*N 
			   		list(c(dN))
			   	})
			   },
			   parms = Parameters,
			   times = c(from = 0, to = 4, by = 0.1), # the time interval over which the model will be simulated.
			   init = c(N = 10), # starting biovolume
			   solver = "lsoda" #lsoda will be called with tolerances of 1e-9. Default tolerances are both 1e-6. Lower is more accurate.
)

fittedparms <- c("r", "alpha") # for assigning fitted parameter values to fittedCRmodel

LowerBound <- c(r = 0.01, alpha = 0.001)
UpperBound <- c(r = 3, alpha = 0.1) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound


controlfit <- function(data){
	
	init(CRmodel) <- c(N = data$N[1]) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, N) # The Y values of the observed data points we are fitting our model to
	
	
	fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
								 debuglevel = 0, fn = ssqOdeModel,
								 method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
								 control = list(trace = T)
	)
	
	r <- coef(fittedCRmodel)[1]
	alpha <- coef(fittedCRmodel)[2]
	well_plate <- data$well_plate[1]
	output <- data.frame(well_plate, r, alpha)
	return(output)
}

# fit models --------------------------------------------------------------

Parameters <- c(r = 1, alpha = 0.01) ## initial parameter guesses

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.01, alpha = 0.001)
UpperBound <- c(r = 5, alpha = 0.5) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound


data <- read_csv("data-raw/nitrate-abundances-processed.csv") %>% 
	select(days, RFU, nitrate_concentration, population, well_plate) %>% 
	rename(N = RFU) %>% 
	filter(population != "COMBO")

nitrate_key <- read_csv("data-raw/nitrate-abundances-processed.csv") %>% 
	select(nitrate_concentration, population, well_plate) %>% 
	distinct(well_plate, .keep_all = TRUE)


data_split <- data %>% 
	split(.$well_plate)

# data %>% 
# 	group_by(well_plate) %>% 
# 	summarise(max_days = max(days)) %>% View

output_CR_all <- data_split %>% 
	map_df(controlfit)

alphas <- left_join(output_CR_all, nitrate_key) %>% 
	filter(population != "COMBO")

alphas %>% 
	ggplot(aes(x = nitrate_concentration, y = alpha)) + geom_point() +
	facet_wrap( ~ population)

write_csv(alphas, "data-processed/alphas-nitrate.csv")



# Plot the fits -----------------------------------------------------------

alphas_split <- alphas %>%
	split(.$well_plate)

df <- alphas_split[[1]]

simulate_growth <- function(df){
	
	plotfittedCRmodel <- CRmodel
	init(plotfittedCRmodel) <- c(N = data$N[1])
	parms(plotfittedCRmodel)[fittedparms] <- c(df$r[1], df$alpha[1])
	
	# set model parameters to fitted values and simulate again
	times(plotfittedCRmodel) <- c(from=0, to=4, by=0.1)
	ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))
	return(ysim)
	
}  

all_simulated_alphas <- alphas_split %>% 
	map_df(simulate_growth, .id = "well_plate")




all_sim <- all_simulated_alphas %>% 
	left_join(., nitrate_key) %>%
	rename(days = time)
	
	all_sim %>% 
	# filter(grepl("B", well_plate)) %>% 
	# filter(well_plate == "B03_1") %>% 
	ggplot(aes(x = days, y = N, group = well_plate, color = factor(nitrate_concentration))) +
	geom_line() +
	geom_point(aes(x = days, y = N, color = factor(nitrate_concentration)), data= filter(data_nitrate)) +
	facet_grid(nitrate_concentration ~ population, scales = "free") +
	scale_color_viridis_d()
	ggsave("figures/r-alpha-logistic.png", width = 30, height = 20)



alphas %>% 
	ggplot(aes(x = nitrate_concentration, y = r, group = population)) + geom_point() +
	geom_line()


mean_r <- alphas %>% 
	filter(nitrate_concentration == 1000) %>% 
	group_by(population) %>% 
	summarise_each(funs(mean, max), r) %>% 
	select(population, mean) %>% 
	rename(r = mean)

max_r <- alphas %>% 
	# filter(nitrate_concentration == "1000") %>% 
	group_by(population) %>% 
	summarise_each(funs(mean, max), r) %>% 
	select(population, max) %>% 
	rename(r = max)

nitrate_r <- left_join(data_nitrate, mean_r)

### new models with fixed r

Parameters <- c(alpha = 0.01) ## initial parameter guesses
fixed_r_growth <- new("odeModel",
			   main = function (time, init, parms) {
			   	with(as.list(c(init, parms)), {
			   		dN <-  (r - alpha*N)*N 
			   		list(c(dN))
			   	})
			   },
			   parms = Parameters,
			   times = c(from = 0, to = 4, by = 0.1), # the time interval over which the model will be simulated.
			   init = c(N = exp(0.781819)), # starting biovolume
			   solver = "lsoda" #lsoda will be called with tolerances of 1e-9. Default tolerances are both 1e-6. Lower is more accurate.
)

fittedparms <- c("alpha") # for assigning fitted parameter values to fittedCRmodel

LowerBound <- c(alpha = 0.001)
UpperBound <- c(alpha = 0.1) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound

data <- data_split_r[[1]]

fit_r <- function(data){
	
	init(CRmodel) <- c(N = data$N[1]) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, N) # The Y values of the observed data points we are fitting our model to
	r <- data$r[[1]]
	
	fittedCRmodel <- fitOdeModel(fixed_r_growth, whichpar = fittedparms, obstime, yobs,
								 debuglevel = 0, fn = ssqOdeModel,
								 method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
								 control = list(trace = T)
	)
	
	# r <- coef(fittedCRmodel)[1]
	alpha <- coef(fittedCRmodel)[1]
	well_plate <- data$well_plate[1]
	output <- data.frame(well_plate, r, alpha)
	return(output)
}

data_split_r <- nitrate_r %>% 
	split(.$well_plate)

# data %>% 
# 	group_by(well_plate) %>% 
# 	summarise(max_days = max(days)) %>% View

output_fixed_r <- data_split_r %>% 
	map_df(fit_r)

output_fixed_r2 <- left_join(output_fixed_r, nitrate_key)


correct_alphas <- output_fixed_r2 %>% 
	select(-r)

correct_rs <- alphas %>% 
	select(-alpha)

all_r_alpha <- left_join(correct_rs, correct_alphas)
write_csv(all_r_alpha, "data-processed/r-alphas-nitrate.csv")


output_fixed_r2 %>% 
	filter(nitrate_concentration < 100) %>% 
	ggplot(aes(x = nitrate_concentration, y = alpha)) + geom_point() + facet_wrap( ~ population, scales = "free")
ggsave("figures/alpha_N.png", height = 15, width = 20)



# hyperbolic function -----------------------------------------------------


# OK so let's fit a hyperbolic function to these.
library(broom)
library(nls.multstart)

fits_many <- output_fixed_r2 %>% 
	group_by(population) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(alpha ~a/(k + nitrate_concentration),
												  data = .x,
												  iter = 500,
												  start_lower = c(a = 0, k = 0.01),
												  start_upper = c(a = 0.5, k = 1),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(a = 0.00000001, k = 0.000001),
												  upper = c(a = 1, k = 25),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))

info <- fits_many %>%
	unnest(fit %>% map(glance))

# get params
params <- fits_many %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 



ks <- output_fixed_r2 %>% 
	group_by(population) %>% 
	do(tidy(nls(alpha ~ a/(k + nitrate_concentration),
				data= .,  start=list(a = 0.2, k = 0.5), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))


predict_hyperbolic <- function(df) {
	
	hyperboliccurve<-function(x){
		alpha <- (df$a[[1]] / (df$k[[1]] +x))
		alpha}
	
	pred <- function(x) {
		y <- hyperboliccurve(x)
	}
	
	x <- seq(0, 1000, by = 0.1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(nitrate_concentration.x = x, 
			   alpha = preds)
}


bs_split <- params %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_n <- bs_split %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_hyperbolic, .id = "population")


output_fixed_r2 %>% 
	# mutate(estimate = mu) %>%
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	filter(population == 29) %>% 
	ggplot(aes(x= nitrate_concentration, y= alpha)) + 
	geom_point() + 
	geom_line(data= filter(all_preds_n, population == 29), aes(x=nitrate_concentration.x, y=alpha), size = 1)+
	facet_wrap( ~ population, scales = "free") + xlim(0, 80) + ylim(0, 0.05)


params %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate) %>% 
	filter(population == 29) %>% 
	mutate(ak = a/k) %>% View

ggsave("figures/alpha_fits.png", height = 15, width = 20)

library(readxl)
library(janitor)
library(plotrix)
treatments <- read_excel("data-raw/ChlamEE_Treatments_JB.xlsx") %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	filter(population != "cc1629")

params %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate) %>% 
	left_join(., treatments) %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), k) %>% 
	ggplot(aes(x = treatment, y = k, color = ancestor_id)) + geom_point()
	

m<- 0.1
rstars_alpha <- params %>% 
	filter(term == "k") %>% 
	left_join(., mean_r) %>% 
	rename(k = estimate) %>% 
	mutate(rstar_solve = k*m/(r-m)) %>% 
	left_join(., treatments)


rstars_alpha %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar_solve), y = rstar_solve, color = ancestor_id), size = 2, data = rstars_alpha, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor")
ggsave("figures/r-star-alpha.png", width = 6, height = 4)




# find dn nt --------------------------------------------------------------


dn <- function(r, alpha, N){
	dn <- (r - alpha*N)
	return(dn)
}

# alphas %>% 
# 	map_df(dn, .id = "well_plate") %>% View


## look at the Ks
Ks <- output_fixed_r2 %>% 
	# all_r_alpha %>% 
	mutate(K = r/alpha) 

Ks %>% 
	ggplot(aes(x = r, y = alpha, color = population)) + geom_point()


all_r_alpha %>% 
	mutate(dn_1 = r - alpha*1) %>% 
	mutate(dn_10 = r - alpha*10) %>% 
	mutate(dn_50 = r - alpha*50) %>% 
	mutate(dn_100 = r - alpha*100) %>% 
	mutate(K = r/alpha) %>% 
	# mutate(dn_1000 = r - alpha*1000) %>% 
	select(population, nitrate_concentration, contains("dn")) %>% 
	gather(key = N, value = dN, contains("dn")) %>% 
	mutate(N  = case_when(N == "dn_10" ~ 10,
						  N == "dn_50" ~ 50,
						  N == "dn_100" ~ 100,
						  N == "dn_1" ~ 1)) %>% 
	ggplot(aes(x = nitrate_concentration, y = dN, color = N)) + geom_point() + facet_wrap( ~ N) +
	geom_hline(yintercept = 0) + ylab("dN/dT/N") + xlab("Nitrate (uM)") 
ggsave("figures/dNdT_figs_nitrate.png", width = 10, height = 8)

Ks %>% 
	ggplot(aes(x = nitrate_concentration, y = K, color = log(alpha))) + geom_point() +
	xlim(0, 30)

#### now lets's plot the fits

## all_r_alpha these are the parameter values we want to simulate

df <- output_fixed_r2 %>% 
	filter(well_plate == "B02_11")


simulate_growth <- function(df){
plotfittedCRmodel <- CRmodel
parms(plotfittedCRmodel)[fittedparms] <- c(df$r[1], df$alpha[1])

# set model parameters to fitted values and simulate again
times(plotfittedCRmodel) <- c(from=0, to=4, by=0.1)
ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))

return(ysim)

}  

output_fixed_r2_split <- output_fixed_r2 %>%
	split(.$well_plate)

all_simulated <- output_fixed_r2_split %>% 
	map_df(simulate_growth, .id = "well_plate")

all_simulated2 <- all_simulated %>% 
	left_join(., nitrate_key)

all_simulated2 %>%
	rename(days = time) %>% 
	filter(grepl("B", well_plate)) %>% 
	ggplot(aes(x = days, y = N, color = factor(nitrate_concentration), group = well_plate)) + geom_line() +
	# geom_point(aes(x = days, y = N, color = factor(nitrate_concentration)), data= data_nitrate) +
	facet_wrap(~ well_plate) +
	scale_color_viridis_d()
ggsave("figures/all-r-alphas-fits.png", width = 30, height = 30)
