predict_hyperbolic <- function(df) {
	
	hyperboliccurve<-function(x){
		alpha <- (df$a[[1]] / (df$k[[1]] +x))
		growth_rate}
	
	pred <- function(x) {
		y <- monodcurve(x)
	}
	
	x <- seq(0, 1000, by = 0.1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(nitrate_concentration.x = x, 
			   alpha = preds)
}