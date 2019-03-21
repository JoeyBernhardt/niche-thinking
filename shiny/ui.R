library(shiny)
library(shinyWidgets)

shinyUI(fluidPage(
  h4("Linking resource competition and Lotka-Volterra models"),
  h5("Patrick Thompson, Joey Bernhardt, Mary O'Connor"),
  h6(""),
  fluidRow(
    column(2,
           wellPanel(
             sliderInput("c11",
                         "c11 - per capita consumption of species 1 on resource 1",
                         min = 0,
                         max = 5,
                         value = 2),
             sliderInput("c12",
                         "c12 - per capita consumption of species 1 on resource 2",
                         min = 0,
                         max = 5,
                         value = 3)
           )),
    column(2,
           wellPanel(
             sliderInput("c21",
                         "c21 - per capita consumption of species 2 on resource 1",
                         min = 0,
                         max = 5,
                         value = 4),
             sliderInput("c22",
                         "c22 - per capita consumption of species 2 on resource 2",
                         min = 0,
                         max = 5,
                         value = 2)
           )),
    column(2,
           wellPanel(
             sliderInput("w11",
                         "w11 - conversion efficiency of species 1 for resource 1",
                         min = 0,
                         max = 5,
                         value = 2),
             sliderInput("w12",
                         "w12 - conversion efficiency of species 1 for resource 2",
                         min = 0,
                         max = 5,
                         value = 4)
           )),
    column(2,
           wellPanel(
             sliderInput("w21",
                         "w21 - conversion efficiency of species 2 for resource 1",
                         min = 0,
                         max = 5,
                         value = 4),
             sliderInput("w22",
                         "w22 - conversion efficiency of species 2 for resource 2",
                         min = 0,
                         max = 5,
                         value = 2)
           )),
    column(2,
           wellPanel(
             sliderInput("S1",
                         "S1 - supply of resource 1",
                         min = 0,
                         max = 50,
                         value = 20),
             sliderInput("S2",
                         "S2 - supply of resource 2",
                         min = 0,
                         max = 50,
                         value = 20)
           ))
    ,
    column(2, 
           wellPanel(
             textOutput("alpha_text"),
             textOutput("alpha_text2"),
             textOutput("co_exist_text1"),
             textOutput("co_exist_text2"),
             textOutput("co_exist_text3")
           ))
  ),
  fluidRow(
    column(3,offset = 0,
           plotOutput("res_comp_plot", width="100%")
    ),
    column(3,
    plotOutput("ZINGI_plot")
    ),
    
    column(3, offset = 0,
           plotOutput("LV_plot", width="100%")
    ),
    column(3, offset = 0,
           plotOutput("coexistence_plot", width="100%")
  )
)
)
)