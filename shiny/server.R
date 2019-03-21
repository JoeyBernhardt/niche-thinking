## This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Patrick L. Thompson
# patrick.thompson@zoology.ubc.ca
# March 2019

library(shiny)
library(deSolve)
library(ggplot2)
library(cowplot)
library(shinyWidgets)

#continuous time resource competition model


server <- function(input, output) {
  model_outputs <- reactive({
    res_comp <- function(t, State, Pars) 
    {
      with(as.list(c(State, Pars)),
           {
             dN1dt <- r1 * N1 * ((w11 * R1 + w12*R2 - T1)/(k1 + w11 * R1 + w12 * R2 - T1))- D * N1 
             dN2dt <- r2 * N2 * ((w21 * R1 + w22*R2 - T2)/(k2 + w21 * R1 + w22 * R2 - T2))- D * N2 
             dR1dt <- D * (S1 - R1) - c11 * N1 - c21 * N2
             dR2dt <- D * (S2 - R2) - c12 * N1 - c22 * N2
             list(c(dN1dt, dN2dt, dR1dt, dR2dt))
           })
    }
    
    #continuous Lotka-Volterra competition model
    LV=function(Time, State, Pars)
    {
      with(as.list(c(State, Pars)),
           {
             dN1 = r1*N1*(1 - a11*N1 - a12*N2)
             dN2 = r2*N2*(1 - a21*N1 - a22*N2)
             return(list(c(dN1, dN2)))
           })
    }
    
    
    D = 0.7 #mortality rate
    c11 = input$c11; c12 = input$c12; c21 = input$c21; c22 = input$c22 ### cij = per capita consumption of comsumer i on resource j
    w11 = input$w11; w12 = input$w12; w21 = input$w21; w22 = input$w22 ## wij = weighting factor that converts availability of resource j into consumer Ni
    k1 = 0.1; k2 = 0.1 #half saturation constant for N resource consumption
    r1 = 0.8; r2 = 0.8 #population growth rates
    T1 = 0.1; T2 = 0.1 #minimum amount of resource required for growth
    
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
    
    #' Supply concentrations of resources
    S1 <- input$S1
    S2 <- input$S2
    
    ZNGI.df <-  data.frame(blue = c(R1.star), red = c(R2.star))
    params_df <- data.frame(c11, c21, c22, c12, R1.star, R2.star, S1 = S1, S2 = S2, y.inter.1, y.inter.2, slope.1, slope.2)
    
    #initial pop size
    N1 <- 0.5
    N2 <- 0.1
    
    #Calculate alphas
    a11 <- (c12 + c11 * (w11/w12)) / (D * (S2 + (w11/w12) * S1 - B1))
    a12 <- (c22 + c21 * (w11/w12)) / (D * (S2 + (w11/w12) * S1 - B1))
    a22 <- (c22 + c21 * (w21/w22)) / (D*(S2 + (w21/w22) * S1 - B2))
    a21 <- (c12 + c11 * (w21/w22)) / (D*(S2 + (w21/w22) * S1 - B2))
    
    alphas <- data.frame(a11 = a11, a12 = a12, a22 = a22, a21 = a21)
    
    #calculate stabilizing potential and fitness ratio
    p <- sqrt((a12*a21)/(a11*a22)) #niche overlap
    stabil_potential <- 1 - p #stabilizing potential
    fit_ratio <- sqrt((a11*a12)/(a22*a21)) #fitness ratio
    coexist <- p < fit_ratio &  fit_ratio < 1/p #coexistence?
    
    coexistence <- data.frame(stabil_potential = stabil_potential, fit_ratio = fit_ratio, coexist = coexist)
    
    parms_rc <- c(r1 = r1, r2 = r2, w11 = w11, w12 = w12, w21 = w21, w22 = w22, T1 = T1, T2 = T2, D = D, c11 = c11, c12 = c12, c21 = c21, c22 = c22)
    initial_rc <- c(N1 = N1, N2 = N2, R1 = 1, R2 = 1)
    
    times = seq(1,200, by = 0.1)
    
    out_rc <- data.frame(ode(y = initial_rc, times = times, func = res_comp, parms = parms_rc, atol = 0))
    
    parms_LV = c(r1=r1, r2=r2, a11=a11, a21=a21, a12=a12, a22=a22)
    initial_LV = c(N1 = N1, N2 = N2)
    
    out_LV = data.frame(ode(y = initial_LV, times = times, func = LV, parms = parms_LV, method = "ode45"))

    list(out_rc = out_rc, out_LV = out_LV, alphas = alphas, coexistence = coexistence, ZNGI.df = ZNGI.df, params_df = params_df)
  })
  
  output$alpha_text <- renderText({ 
    hold <- model_outputs()
    alphas <- hold$alphas
    paste("a11 = ", round(alphas$a11,3),", a12 = ", round(alphas$a12,3))
  })
  
  output$alpha_text2 <- renderText({ 
    hold <- model_outputs()
    alphas <- hold$alphas
    paste("a22 = ", round(alphas$a22,3),", a21 = ", round(alphas$a21,3))
  })
  
  output$co_exist_text1 <- renderText({ 
    hold <- model_outputs()
    coexistence <- hold$coexistence
    paste("stabilizing potential = ", round(coexistence$stabil_potential,2))
  })
  
  output$co_exist_text2 <- renderText({ 
    hold <- model_outputs()
    coexistence <- hold$coexistence
    paste("fitness ratio = ", round(coexistence$fit_ratio,2))
  })
  
  output$co_exist_text3 <- renderText({ 
    hold <- model_outputs()
    coexistence <- hold$coexistence
    paste("coexistence = ", coexistence$coexist)
  })
  
  
  output$res_comp_plot <- renderPlot({
    hold <- model_outputs()
    out_rc <- hold$out_rc
    out_LV <- hold$out_LV
    plot(out_rc$time, out_rc$N1, type = "n", xlab = "time", ylab = "population size", ylim = c(0,max(out_rc[,2:3])), main = "resource competition model")
    abline(h = out_rc$N1[nrow(out_rc)], lty = 2, col = 2)
    abline(h = out_rc$N2[nrow(out_rc)], lty = 2, col = "dodgerblue")
    lines(out_rc$time, out_rc$N1, col = 2, lwd = 2)
    lines(out_rc$time, out_rc$N2, col = "dodgerblue", lwd = 2)
    legend("bottomright", col = c(2,"dodgerblue"), lty = 1, lwd = 2, legend = c("species 1", "species 2"), bty = "n")
  })
  
  output$LV_plot <- renderPlot({
    hold <- model_outputs()
    out_rc <- hold$out_rc
    out_LV <- hold$out_LV
    plot(out_LV$time, out_LV$N1, type = "n", xlab = "time", ylab = "population size", ylim = c(0,max(out_LV[,2:3])), main = "Lotka-Volterra competition model")
    abline(h = out_LV$N1[nrow(out_LV)], lty = 2, col = 2)
    abline(h = out_LV$N2[nrow(out_LV)], lty = 2, col = "dodgerblue")
    lines(out_LV$time, out_LV$N1, col = 2, lwd = 2)
    lines(out_LV$time, out_LV$N2, col = "dodgerblue", lwd = 2)
    legend("bottomright", col = c(2,"dodgerblue"), lty = 1, lwd = 2, legend = c("species 1", "species 2"), bty = "n")
  })
  
  output$coexistence_plot <- renderPlot({
    hold <- model_outputs()
    coexistence <- hold$coexistence
    alphas <- hold$alphas 
    
    a11 <- alphas$a11;  a12 <- alphas$a12 ;  a22 <- alphas$a21;  a11 <- alphas$a21
    
    x.1 = seq(1, 0, by=-0.001)
    x.2 = seq(0, -1, by=-0.001)
    y1.1 = 1/((1-x.1))
    y2.1 = 1-x.1
    y1.2 = 1/((1-x.2))
    y2.2 = 1-x.2
    
    cf.df.1 = data.frame(my.x = c(x.1, x.1), 
                         my.y = c(y1.1, y2.1), 
                         whichy = c(rep("y1", times = length(seq(1, 0, by=-0.001))), 
                                    rep("y2", times = length(seq(1, 0, by=-0.001)))))
    rib.dims.1 = data.frame(min.dim = y1.1, max.dim = y2.1, x.dim = x.1)
    
    cf.df.2 = data.frame(my.x = c(x.2, x.2), 
                         my.y = c(y1.2, y2.2), 
                         whichy = c(rep("y1", times = length(seq(0, -1, by=-0.001))), 
                                    rep("y2", times = length(seq(0, -1, by=-0.001)))))
    rib.dims.2 = data.frame(min.dim = y1.2, max.dim = y2.2, x.dim = x.2)
    
    
    ggplot(coexistence, aes(x = 1, y = 1))+
      geom_ribbon(data = rib.dims.1, 
                  aes(x = x.dim, 
                      ymin = min.dim, 
                      ymax = max.dim), 
                  fill = "grey", 
                  alpha = 0.3) +
      geom_line(data = cf.df.1, 
                aes(x = my.x, 
                    y = my.y,
                    linetype = whichy,
                    group = whichy), 
                col = "black") +     
      geom_ribbon(data = rib.dims.2, 
                  aes(x = x.dim, 
                      ymin = min.dim, 
                      ymax = max.dim), 
                  fill = "darkgrey", 
                  alpha = 0.6) +
      geom_line(data = cf.df.2, 
                aes(x = my.x, 
                    y = my.y,
                    linetype = whichy,
                    group = whichy), 
                col = "black") + 
      geom_abline(intercept = 1, 
                  slope = 0, 
                  lty = 2, 
                  col = "darkgrey") + 
      geom_vline(xintercept = 0, 
                 lty = 2, 
                 col = "darkgrey") +
      panel_border(colour = "black") +
      scale_linetype_manual(values = c("y2" = "solid", 
                                       "y1" = "dotted")) +
      xlab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) + 
      ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
      theme(legend.position = "none", 
            axis.title.y = element_text(angle = 90), 
            axis.title = element_text(size = 20))+
      coord_cartesian(ylim = c(0,2))+
      geom_text(aes(x = 0.0,
                    y = 0.25),
                label = "Species 1 wins",
                size = 6.0) +  
      geom_text(aes(x = 0.0,
                    y = 1.75),
                label = "Species 2 wins",
                size = 6.0) +  
      geom_text(aes(x = 0.7,
                    y = 1),
                label = "Coexistence",
                size = 6.0) +  
      geom_text(aes(x = -0.7,
                    y = 1),
                label = "Priority effect",
                size = 6.0)+
      geom_point(data = coexistence, aes(x = stabil_potential, y = fit_ratio), size = 3, fill = "red", pch = 21)
    })
  
  output$ZINGI_plot <- renderPlot({
    hold <- model_outputs()
    ZNGI.df <- hold$ZNGI.df
    params_df <- hold$params_df
    
    ggplot(ZNGI.df, aes(x=blue, y=red)) +   
      geom_point(x = params_df$S1/(2*max(params_df$S1, params_df$S2)), y = params_df$S2/(2*max(params_df$S1, params_df$S2)), size = 3) +
      geom_abline(data = params_df, aes(intercept = y.inter.1, slope = slope.1), col ='#d1495b', size = 0.5) + ## red species' ZNGI
      geom_abline(data = params_df, aes(intercept = y.inter.2, slope = slope.2), col ='#30638e', size = 0.5) + ## blue species' ZNGI
      coord_cartesian(expand = 0, ylim=c(0, 0.6), xlim = c(0, 0.6)) +
      
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
     })
}