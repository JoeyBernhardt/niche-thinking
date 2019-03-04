
#' ---
#' title: "Playing around with ZNGIs"
#' author: "Joey"
#' ---


#+
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(cache = TRUE)

#+ knitr::opts_chunk$set(message = FALSE)


library(tidyverse)
library(cowplot)

#' Setup Tilman's consumer resource model parameters 

### Here species 1 is RED species 2 is BLUE

D = 0.7
c11 = 2; c12 = 4; c21 = 4; c22 = 2 ### cij = per capita consumption of comsumer i on resource j
w11 = 2; w12 = 4; w21 = 4; w22 = 2 ## wij = weighting factor that converts availability of resource j into consumer Ni
k1 = 0.4; k2 = 0.4
r1 = 1; r2 = 1
T1 = 0.1; T2 = 0.1

species1_consumption <- c12/c11
species2_consumption <- c22/c21

#' Calculate slope and intercept for the ZNGI's (get this from equation S3 in Ke and Letten)
slope.1 <-  - (w11 / w12)
slope.2 <-  - (w21 / w22)
y.inter.1 <-  (D * (k1 - T1) + r1 * T1) / (w12 * (r1 - D))
y.inter.2 <-  (w21 / w22) * (D * (k2 - T2) + r2 * T2) / (w21 * (r2 - D))


#' Calculate equilibrium resource values (when the ZNGIs cross each other, i.e. such that Ke and Letten equations S3.1 = S3.2)
B1 <-  y.inter.1
B2 <-  y.inter.2
Lamda.1 <-  (w11 / w12)
Lamda.2 <- (w21 / w22)
R1.star <-  (B1 - B2) / (Lamda.1 - Lamda.2)
R2.star <-  (B2 * Lamda.1 - B1 * Lamda.2) / (Lamda.1 - Lamda.2)  


#' Supply points
S1 <-  0.248; S2 <-  0.245

# plot settings
mylims.x <-  0.55
mylims.y <-  0.55


#' get all the parameters together
ZNGI.df <-  data.frame(blue = c(R1.star), red = c(R2.star))
params_df <- data.frame(c11, c21, c22, c12, R1.star, R2.star)


#' plot it!
ggplot(ZNGI.df, aes(x=blue, y=red)) +   
	geom_point(x = S1, y = S2, size = 3) +
	geom_abline(intercept = y.inter.1, slope = slope.1, col ='#d1495b', size = 0.5) + ## red species' ZNGI
	geom_abline(intercept = y.inter.2, slope = slope.2, col ='#30638e', size = 0.5) + ## blue species' ZNGI
	coord_cartesian(expand = 0, ylim=c(0, mylims.y), xlim = c(0, mylims.x)) +
	
	geom_segment(data = params_df, 
				 aes(x = R1.star, y = R2.star, xend = R1.star-c11/100, yend = R2.star-c12/100),
				 size = 0.5, 
				 col='#d1495b', 
				 arrow = arrow(type = "closed", length = unit(0.05, "inches"))) + 
	geom_segment(data = params_df, 
				 aes(x = R1.star, y = R2.star, xend = R1.star-c21/100, yend = R2.star-c22/100), 
				 size = 0.5, col='#30638e', 
				 arrow = arrow(type = "closed", length = unit(0.05, "inches"))) + 
	geom_segment(data = params_df, 
				 aes(x = R1.star, y = R2.star, xend = R1.star+c11, yend = R2.star+c12), 
				 size = 0.5, 
				 col='#d1495b', 
				 linetype = 2) + 
	geom_segment(data = params_df, 
				 aes(x = R1.star, y = R2.star, xend = R1.star+c21, yend = R2.star+c22), 
				 size = 0.5, 
				 col='#30638e', 
				 linetype = 2) + 
	xlab(expression(R[1])) + 
	ylab(expression(R[2])) +
	
	theme(legend.position = "none", 
		  plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		  axis.text = element_text(size=13),
		  axis.title=element_text(size=20)) +
	panel_border(colour = "black")
