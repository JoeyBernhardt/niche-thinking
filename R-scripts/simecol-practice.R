

### experimenting with the chemostat data in the simecol package


library(simecol)


data("chemostat")
solver(chemostat)


solver(chemostat) <- function(y, times, func, parms, ...) {
	 ode(y, times, func, parms, method = "adams", ...)
}


cs1 <- cs2 <- chemostat


str(chemostat)


obstime <- seq(0, 20, 2)


yobs <- data.frame(
	X = c(10, 26, 120, 197, 354, 577, 628, 661, 654, 608, 642),
	S = c(9.6, 10.2, 9.5, 8.2, 6.4, 4.9, 4.2, 3.8, 2.5, 3.8, 3.9)
	)

yobs <- select(data, N, R) %>% 
	rename(X = N,
		   S = R)
obstime <- data$days
times(cs1) <- obstime

yobs %>% 
	ggplot(aes(x = obstime, y = X)) + geom_point() +
	geom_point(aes(x = obstime, y = yobs$S), color = "red")


res <- fitOdeModel(cs1, obstime = obstime, yobs=yobs)
warnings()


whichpar <- c("vm", "km", "Y")
lower <- c(vm=0, km=0, Y=0)
upper <- c(vm=100, km=500, Y=200)
parms(cs1)[whichpar] <- c(vm=5, km=10, Y=100)
res <- fitOdeModel(cs1, whichpar = whichpar,
					  lower = lower, upper=upper,
					  obstime = obstime, yobs = yobs, method = "PORT",
					  control=list(trace = FALSE))
res


### try with my data

data <- read_csv("data-raw/population4-nitrate.csv") %>% 
	select(days, RFU, nitrate_concentration, population, well_plate) %>% 
	rename(starting_nitrate = nitrate_concentration) %>% 
	mutate(nitrate = starting_nitrate - RFU/20) %>% 
	rename(N = RFU,
		   R = nitrate) %>% 
	filter(starting_nitrate == 1000, well_plate == "F03_29") 

yobs <- select(data, N, R) %>% 
	rename(X = N,
		   S = R)
obstime <- data$days
times(cs1) <- obstime

whichpar <- c("vm", "km")
pars <- c("vm", "km", "Y", "S0", "D")
lower <- c(vm=0, km=0, Y=0)
upper <- c(vm=100, km=500, Y=200)
# parms(cs1)[whichpar] <- c(vm=1, km=10, Y=1)
parms(cs1)[pars] <- c(vm=1, km=10, Y=1, S0 = 999.25, D = 0)

init(cs1) <- c(X = data$N[1], S = data$R[1])

res <- fitOdeModel(cs1, whichpar = whichpar,
				   lower = lower, upper=upper,
				   obstime = obstime, yobs = yobs, method = "PORT",
				   control=list(trace = FALSE))
res

yobs <- select(data, N, R) %>% 
	rename(X = N,
		   S = R)
obstime <- data$days
times(cs1) <- obstime

yobs %>% 
	ggplot(aes(x = obstime, y = X)) + geom_point() +
	geom_point(aes(x = obstime, y = yobs$S), color = "red")