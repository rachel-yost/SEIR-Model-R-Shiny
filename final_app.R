#Required libraries
library(deSolve)
library(tidyr)
library(reshape)
library(magrittr)
library(plyr)
library(dplyr)
library(ggplot2)
library(shiny)
library(shinythemes)
library(shinyjs)
library(plotly)

m <- 3  #number of age classes
mu <- c(0,0,0.00) #death rate in each age group, we are leaving the death rate at zero
S0 <- c(.99,0.99,0.99) # initial value for number of susceptible per age group
E0 <- c(0.005,0.005,0.005) # initial value for number of exposed per age group
I0 <- c(0.004,0.004,0.004) # initial value for number of infectious per age group
R0 <- c(0.00, 0.00,.00) # initial value for number of recovered per age group
num_days <- 365 # number of days to simulate

# Define UI for application
ui <- fluidPage(
  useShinyjs(), # initialize package that allows the action button to start already pressed
  theme=shinytheme("slate"),
  titlePanel("Interactive Age-Structured SEIR model"),
  br(),
  uiOutput("text1"), #paragraph of text, using uiOutput so unicode symbols can be displayed
  br(),
  # Adding buttons to choose parameter values
  sidebarPanel(
    numericInput("beta1input", HTML("&beta;1"), 1,
                 0, 20, .01),
    numericInput("beta2input", HTML("&beta;2"), 1,
                 0, 20, .01),
    numericInput("beta3input", HTML("&beta;3"), 1,
                 0, 20, .01),
    numericInput("sigmainput", HTML("&sigma;"), 0.179,
                 0, 20, .01),
    numericInput("gammainput", HTML("&gamma;"), .5,
                 0, 20, .01),
    selectInput("AgeInput", "Age Groups",c("All",
                                           "0-19 years",
                                           "19-49 years",
                                           "50-75 years")),
    actionButton("gobutton", "Get SEIR Model") #when user clicks, plots and stats are displayed
  ),
  mainPanel(
    tabsetPanel(tabPanel("Plot", plotlyOutput("plotly")),
                tabPanel("Statistics", tableOutput("resultsTable")),
                tabPanel("References",verbatimTextOutput("References"))
    )
  )
)

server <- function(input, output) {
  SEIR_model <- reactive({
    #creating matrix of beta coefficients
    b11 <- input$beta1input
    b22 <- input$beta2input
    b33 <- input$beta3input
    #the between-group beta is determined by the smaller beta value
    b12 <- ifelse(input$beta1input < input$beta2input, input$beta1input, input$beta2input)
    b13 <- ifelse(input$beta3input < input$beta1input, input$beta3input, input$beta1input)
    b21 <- ifelse(input$beta1input < input$beta2input, input$beta1input, input$beta2input)
    b23 <- ifelse(input$beta2input < input$beta3input, input$beta2input, input$beta3input)
    b32 <- ifelse(input$beta2input < input$beta3input, input$beta2input, input$beta3input)
    b31 <- ifelse(input$beta3input < input$beta1input, input$beta3input, input$beta1input)
    
    # parameters
    beta <- matrix(c(b11,b12,b13,b21,b22,b23,b31,b32,b33), nrow=3, ncol=3) 
    gamma <- input$gammainput # recovery rate
    sigma <- input$sigmainput # rate at which individuals move from the exposed to the infectious classes
    
    #creating list of parameters and vector of initial values for the differential equations
    parms <- list(beta=beta, mu=mu, sigma=sigma, gamma=gamma)
    initial <- c(S0, E0, I0, R0)
    t_range <- seq(from= 0, to=num_days, by=1) # vector form 0 to 365 with step of 1
    
    #function to use in the ODE solver
    diff_eqs <- function(times, Y, parms){
      #' Computes the derivatives for the differential equations for each compartment for each age group
      #' for use in the lsoda function
      #' @param times numeric, Current time point
      #' @param Y current variables in the ODE system
      #' @param parms list of parameters used in the equations
      #' @return returns a list of the derivatives for each compartment for each age group
  
      dY <- numeric(length(Y))
      with(parms,{
        for (i in 1:m) { #for each age group, 
          dY[i] <-   -beta[,i]%*%Y[2*m + seq(1:m)] * Y[i] - mu[i] * Y[i] # S_i equation
          dY[m+i] <-  beta[,i] %*% Y[2*m + seq(1:m)] *Y[i] - mu[i] * Y[m+i] - sigma * Y[m+i] #E_i equation
          dY[2*m+i] <- sigma * Y[m+i] - gamma * Y[2*m + i] - mu[i] * Y[2*m+i] #I_i equation
          dY[3*m+i] <- gamma * Y[2*m+i] - mu[i] * Y[3*m + i] #R_i equation
        }
        return(list(dY))
      })
    }
    result = lsoda(initial, t_range, diff_eqs, parms) #solves the differential equations
    #labeling of the output from ODE solver
    compartment <- c("S", "S", "S","E", "E", "E", "I", "I", "I", "R", "R", "R")
    age_group <- c("1", "2", "3","1", "2", "3", "1", "2", "3", "1", "2", "3")
    df <- data.frame(time = result[, 1],
                     compartment= rep(compartment, each = nrow(result)),
                     age_group = rep(age_group, each =  nrow(result)),
                     value = c(result[, -1]))
    #plotting the data
    df$compartment <- factor(df$compartment, levels = c("S","E","I","R"))
    df$age_group <- factor(df$age_group)
    df <- df %>% mutate(compartment = case_when(compartment == "S" ~ "Susceptible",
                                                compartment == "E" ~ "Exposed",
                                                compartment == "I" ~ "Infectious",
                                                compartment == "R" ~ "Recovered"),
                        age_group = case_when(age_group == "1" ~ "0-19 years",
                                              age_group == "2" ~ "19-49 years",
                                              age_group == "3" ~ "50-75 years"))
    if(input$AgeInput == "All") { #if the user wants all age groups
      df
    }
    else { #if the user has filtered the data, only display the df with the selected age group
      df %>% filter(age_group == input$AgeInput)
    }
  }) %>%
    bindEvent(input$gobutton) #make it so this doesn't reload until the action button is pressed
  
  output$plotly <- renderPlotly({
    if(length(unique(SEIR_model()$age_group)) == 3) { #if the user selected all ages:
      ggplot(data = SEIR_model()) + geom_line(aes(x = time, y = value, color = age_group)) +
        theme_minimal()+
        facet_wrap(~factor(compartment, levels = c("Susceptible", "Exposed", "Infectious", "Recovered")), ncol=1, scales =  "free_y") +
        xlab("Time (days)") + ylab("Proportion of Population") +
        guides(color = guide_legend(title = "Age Group"))
    }
    else{ #if the user selected a specific age group
      plot_ly(SEIR_model(), x = ~time, y = ~value, color = ~compartment, type = 'scatter', mode = 'lines+markers') %>%
        layout(
          xaxis = list(
            rangeslider = list(type = "Time (days)")
          ),
          yaxis = list(title = "Proportion of Population")
        )
    }
  })
  output$resultsTable <- renderTable({
    SEIR_model() %>% group_by(age_group, compartment) %>%
      filter(compartment != "Susceptible") %>% #remove susceptible compartment bc the max will always be .99
      summarize(max = max(value)) %>%
      rename("Age Group" = age_group, "Compartment" = compartment,
                     "Max Proportion" = max)
  })
  output$text1 <- renderUI({
    HTML("An SEIR model cateogrizes the host of a disease into one of four compartments.
    Susceptible (unexposed to the pathogen), Exposed (exposed to the pathogen),
    Infected, or Recovered. The parameters involved in forming this model are \u03b2 
     (contact rate), \u03c3 (latency), and \u03b3 (recovery rate). In an
    age-structured SEIR model, the population is
    divided into age categories that have their own dynamics. </p>
    In this SEIR model builder, \u03b21 refers to the contact rate between people under 19 with others
    in that age group, \u03b22 refers to the contact rate between people from 19-49 with others in that age group, and
    \u03b23 is the contact rate between people 50-75 with others in that age group. For our model, the \u03b2 coefficient
    for the group with the lower in-group contact rate determines the between-group contact rate. We based our model off
    of program 3.4 in \"Modeling Infectious Disease in humans and animals\" by Keeling and Rohani.")
  })
  # Auto clicks the play button upon first loading of site
  start_pressed <- observe({
    shinyjs :: click("gobutton")
    start_pressed$destroy() #destroys object so it only auto clicks once
  })
  output$References <-renderText("1. SEIR - SEIRS model. (n.d.). SEIR and SEIRS models. Retrieved December 7, 2023, from <https://docs.idmod.org/projects/emod-hiv/en/latest/model-seir.html>
2. Wickham, H., Vaughan. D. and Girlich, M. (2023). _tidyr: Tidy Messy Data_. R package version 1.3.0, <https://CRAN.R-project.org/package=tidyr>
3. Wickham, H. (2007). Reshaping data with the reshape package. Journal of Statistical Software, 21(12).
4. Bache, S., Wickham, H. (2022). _magrittr: A Forward-Pipe Operator for R_. R package version 2.0.3, <https://CRAN.R-project.org/package=magrittr>.
5. Wickham, H. (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL https://www.jstatsoft.org/v40/i01/.
6. Wickham, H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.
7. Chang, W., Cheng, J., Allaire, J. (2023). _shiny: Web Application Framework for R_. R package version 1.7.5.1, <https://CRAN.R-project.org/package=shiny>.
8. Chang, W. (2021). _shinythemes: Themes for Shiny_. R package version 1.2.0, <https://CRAN.R-project.org/package=shinythemes>.
9. Wickham H., François R., Henry L. (2023). _dplyr: A Grammar of Data Manipulation_. R package version 1.1.3, <https://CRAN.R-project.org/package=dplyr>.
10. Soetaert, K., Petzoldt, T. and Setzer, R. W. (2016). deSolve: General solvers for initial value problems of ordinary differential equations (oDE), partial differential equations (pDE), differential algebraic equations (dAE) and delay differential equations (dDE) [Internet]. Retrieved December 7, 2023, from <https://rdrr.io/cran/deSolve/man/deSolve.html>.
11. Attali, D. (2021). _shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds_. R package version 2.1.0, <https://CRAN.R-project.org/package=shinyjs>.
12. Sievert, C. (2020). Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC. <https://plotly-r.com>
13. Peng, L., Yang, W., Zhang, D., Zhuge, C. and Hong, L. (2020). Epidemic analysis of COVID-19 in China by dynamical modeling. MedRxiv Epidemiol. doi: 10.1101/2020.02.16.20023465.
14. Godio, A., Pace, F. and Vergnano, A. (2020). SEIR Modeling of the Italian Epidemic of SARS-CoV-2 Using Computational Swarm Intelligence. Int J Environ Res Public Health. May 18;17(10):3535. doi: 10.3390/ijerph17103535. PMID: 32443640; PMCID: PMC7277829.
15. U.S. Census Bureau. (n.d.). Annual Estimates of the Resident Population for Selected Age Groups by Sex for the United States: April 1, 2020 to July 1, 2022 (NC-EST2022-AGESEX) Retrieved December 7, 2023, from <https://www.census.gov/data/tables/time-series/demo/popest/2020s-national-detail.html>
16. Bjørnstad, O.N., Shea, K., Krzywinski, M. et al. (2020) The SEIRS model for infectious disease dynamics. Nat Methods 17, 557–558. <https://doi.org/10.1038/s41592-020-0856-2>
17. Keeling, M., Rohani, P. (2008). Modeling Infectious Diseases in humans and animals(pp.82-88). Princeton University Press.
18. Gleeson, P. (2021). How to model an epidemic with R. freeCodeCamp.org. <https://www.freecodecamp.org/news/how-to-model-an-epidemic-with-r/>
19. Stocks, T. (2018) Program 3.4 from Keeling and Rohani. epirecipes. <http://epirecip.es/epicookbook/chapters/kr08/3_4/r_desolve>")
}
# Run the application
shinyApp(ui = ui, server = server)