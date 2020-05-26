#' Create interactive Shiny app to view Observed/Expected output over time by geography and by age
#'
#' \code{dashboardPlotOe} Creates an interactive Shiny plot to explore Observed/Expected 
#'   results generated in the function \code{\link{excessCases}}. Drop down menus
#'   allow for viewing Observed/Expected by different age groups and geographies.
#'
#' @param excess_output Provide the object returned by the function
#'   \code{\link{excessCases}}
#'   
#' @param datevar Provide the datevar used in the function
#'   \code{\link{excessCases}}
#'   
#' @param agevar Provide the agevar used in the function
#'   \code{\link{excessCases}}
#'   
#' @param statevar Provide the statevar geographic variable used in the function
#'   \code{\link{excessCases}} 
#'   
#' @param outcome Provide one of the outcomes / syndromes of interest
#'   \code{\link{excessCases}}  
#'   
#' @param yaxis Provide the variable intended for the y-axis (either agevar or statevar)
#'   \code{\link{excessCases}}  
#'   
#' @param facet Provide the variable intended for the facet groupings (either agevar or statevar)
#'   \code{\link{excessCases}} 
#'   
#' @return Returns an object representing the shinyapp. Depending on the
#'   environment (RStudio vs. console) the app may then be passed to
#'   \code{print()}, which will start the server. Refer to
#'   \code{link[shiny]{shinyApp}} for more details.
#'
#' @export

dashboardPlotOe <- function(excess_output,
                            datevar,  
                            agevar,    
                            statevar,  
                            outcome, 
                            yaxis, 
                            facet) { 
  
  library(dplyr)
  
  #-------------------------------
  #------- Data 
  #-------------------------------
  df_y <- 
    excess_output  %>% 
    pluck(outcome) %>% 
    map(~(
      transpose(.x) %>% 
        map(~(.x %>% bind_rows())) %>% .[c(datevar, "y")] %>% 
        map(~(.x %>% gather(statevar, y))) %>% 
        bind_cols()  %>% 
        set_names(c(statevar, datevar, "rm", "y")) %>% 
        dplyr::select(-rm)
    )) %>% 
    bind_rows(.id = agevar) 
  
  df_pred <-
    excess_output  %>%
    pluck(outcome) %>%
    map(~(
      transpose(.x) %>%
        map(~(.x %>% bind_rows())) %>% .[c(datevar, "pred")] %>%
        map(~(.x %>% gather(statevar, pred))) %>%
        bind_cols()  %>%
        set_names(c(statevar, datevar, "rm", "pred")) %>%
        dplyr::select(-rm)
    )) %>%
    bind_rows(.id = agevar)
  
  df_oe <-
    df_y %>%
    left_join(df_pred, by=c(agevar, statevar, datevar)) %>%
    mutate(oe = y/pred,
           year = year(get(datevar)),
           week = week(get(datevar)),
           oe_fac = cut(oe, 
                        breaks=c(-Inf,
                                 0.5, 1.0,
                                 1.25, 1.5, 1.75, 2.0,
                                 2.5, 3.0,3.5, 4.0,
                                 6.0, 8.0, 10.0,
                                 Inf),
                        labels=c("0.5","1.0",
                                 "1.25", "1.5", "1.75", "2.0",
                                 "2.5", "3.0","3.5", "4.0",
                                 "6.0", "8.0","10.0", ">10.0")),
           oe_fac_rev = factor(oe_fac, levels = rev(levels(oe_fac)))) 
  
  dates <- as.Date(unique(df_oe[[datevar]]))
  states <- unique(df_oe[[statevar]])
  age_groups <- unique(df_oe[[agevar]])
  last.date <- max(dates)
  last.date.format <- format(last.date, "%b %d, %Y")
  
  #-----------------------------------
  #------- UI 
  #-----------------------------------
  ui <- fluidPage(
    
    shiny::titlePanel(paste0('Data through ', last.date.format)),
    
    shiny::sidebarLayout(
      
      shiny::sidebarPanel(
        
        shiny::selectInput(input = "set.states",
                           label = "State:",
                           choice = states,
                           selected = c("NY"),
                           multiple = T),
        
        shiny::selectInput(input = "set.ages",
                           label = "Age group:",
                           choice = age_groups,
                           selected = c("18 and under","19-64","65 and over"),
                           multiple = T),
        
        # note: weekly only
        shiny::sliderInput(
          input = 'display.dates',
          label = 'Earliest date to display',
          min = min(dates),
          max = dates[length(dates)-2],
          step = 7,    
          value = dates[length(dates) - round(length(dates)/5)]
        )),
      
      shiny::mainPanel(shiny::plotOutput("plot")) 
      
    ))
  
  #-----------------------------------
  #-------- Server
  #-----------------------------------
  server <- function(input, output){
    library(ggplot2)
    
    dates_states_ages_sel <- reactive({
      req(input$display.dates, input$set.states, input$set.ages)
      df_oe %>% 
        filter(get(datevar) >= input$display.dates & get(statevar) %in% c(input$set.states) & get(agevar) %in% c(input$set.ages))
    })
    
    output$plot = renderPlot({
      
      ggplot(data = dates_states_ages_sel(), 
             aes(x = factor(get(datevar)), 
                 y = get(yaxis))) +
        geom_raster(aes(fill = oe_fac_rev), interpolate = F) +
        scale_fill_manual(values = c(">10.0"="#5c0900", "10.0"="#850d00",
                                     "8.0"="#a31000", "6.0"="#c21300",
                                     "4.0"="#eb1800", "3.5"="#ff3e29",
                                     "3.0"="#ff7014", "2.5"="#ff9049",
                                     "2.0"="#ffaf35", "1.75"="#ffd230",
                                     "1.5"="#a3a3cc", "1.25"="#b7b7d7",
                                     "1.0"="#cacae2", "0.5"="#dbdbeb")) +
        xlab("Time") +
        labs(fill = "O/E Ratio") +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        facet_grid(get(facet) ~ . ) +
        theme_bw() +
        theme(axis.title.y = element_blank(),
              axis.text.x = element_text(size = 7, vjust = 1, hjust = 0, angle = 90))
    })
    
  }
  
  shiny::shinyApp(ui, server)
  
}
