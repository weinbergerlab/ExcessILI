#' Create interactive Shiny app to explore results
#'
#' \code{dashboardPlot} Creates an interactive Shiny plot to explore results
#'   generated in the function \code{\link{excessCases}}. Drop down menus
#'   allow for viewing different syndromes, age groups, and geographies, and
#'   for looking at plots of raw counts, proportions, or Observed/Expected.
#'
#' @param all.glm.res Provide the object returned by the function
#'   \code{\link{excessCases}}
#'
#' @return Returns an object representing the shinyapp. Depending on the
#'   environment (RStudio vs. console) the app may then be passed to
#'   \code{print()}, which will start the server. Refer to
#'   \code{link[shiny]{shinyApp}} for more details.
#'
#' @export
dashboardPlot <- function(all.glm.res){ 

  ds               <- all.glm.res
  counties.to.test <- names(ds[[1]][[1]])
  ages.to.test     <- names(ds[[1]])
  age.labels       <- ages.to.test
  syndromes        <- names(ds)
  dates            <- as.Date(ds[[1]][[1]][[1]][['date']])
  n.times          <- length(dates)
  last.date.format <- max(dates)
  last.date.format <- format(last.date.format, "%b %d, %Y")

  server <- function(input, output){
    output$countyPlot = renderPlot({

      # Essentially a partialized sapply, one-off for the following few lines
      plucker <-
        function(var) function(x) sapply(x, "[[", var, simplify='array')
      
      data_to_pluck <- ds[[input$set.syndrome]]

      ili2.resid    <- sapply(data_to_pluck, plucker("resid1"), simplify='array')
      ili2.pred     <- sapply(data_to_pluck, plucker("pred"), simplify='array')
      ili2.pred.lcl <- sapply(data_to_pluck, plucker("lpi"), simplify='array')
      ili2.pred.ucl <- sapply(data_to_pluck, plucker("upi"), simplify='array')
      obs.ili       <- sapply(data_to_pluck, plucker("y"), simplify='array')
      denom         <- sapply(data_to_pluck, plucker("denom"), simplify='array')
      
      dimnames(ili2.resid)[[2]]    <- counties.to.test
      dimnames(ili2.pred)[[2]]     <- counties.to.test
      dimnames(ili2.pred.lcl)[[2]] <- counties.to.test
      dimnames(ili2.pred.ucl)[[2]] <- counties.to.test
      dimnames(obs.ili)[[2]]       <- counties.to.test
      dimnames(denom)[[2]]         <- counties.to.test

      plot.min <- which(input$display.dates == dates)
      date.idxs <- plot.min:n.times
      par(mfrow = c(2, 3), mar = c(3, 2, 1, 1))
      
      if(input$arrange.plots=='Age'){
        i.select<-input$set.ages
        j.select<-input$set.borough
      } else {
        i.select<-input$set.ages
        j.select<-counties.to.test
      }

      for(i in i.select){
        for( j in j.select){
            if (input$set.prop == "Counts") {

              y <- obs.ili[date.idxs, j, i]

              pred     <- ili2.pred    [date.idxs, j, i]
              pred.lcl <- ili2.pred.lcl[date.idxs, j, i]
              pred.ucl <- ili2.pred.ucl[date.idxs, j, i]

              if (identical(input$set.axis, F)) {
                maxval <-
                  max(c(ili2.pred.lcl[date.idxs, j, i],
                        ili2.pred.ucl[date.idxs, j, i],
                        ili2.pred    [date.idxs, j, i],
                        obs.ili      [date.idxs, j, i]),
                      na.rm = T)

                y.range <- c(0, maxval)

              } else {
                maxval <-
                  max(c(ili2.pred.lcl[date.idxs, j, ],
                        ili2.pred.ucl[date.idxs, j, ],
                        ili2.pred    [date.idxs, j, ],
                        obs.ili      [date.idxs, j, ]),
                      na.rm = T)

                y.range <- c(0, maxval)

              }
            } else if (input$set.prop == "Proportion") {
              common.denom <- denom[date.idxs, j, i]

              y        <- obs.ili      [date.idxs, j, i] / common.denom
              pred     <- ili2.pred    [date.idxs, j, i] / common.denom
              pred.lcl <- ili2.pred.lcl[date.idxs, j, i] / common.denom
              pred.ucl <- ili2.pred.ucl[date.idxs, j, i] / common.denom
              
              if (input$set.axis == F) {
                y.range <- c(0, max(y, na.rm = T))
              } else {
                y.range <- c(0, max(plot.prop[date.idxs, j, ], na.rm = T))
              }
            } else {
              y        <- obs.ili[date.idxs, j, i]/ili2.pred    [date.idxs, j, i]
              pred     <- obs.ili[date.idxs, j, i]/ili2.pred    [date.idxs, j, i]
              pred.lcl <- obs.ili[date.idxs, j, i]/ili2.pred.lcl[date.idxs, j, i]
              pred.ucl <- obs.ili[date.idxs, j, i]/ili2.pred.ucl[date.idxs, j, i]

              if (input$set.axis == F) {
                y.range <- range(y, na.rm = T)
                y.range[is.infinite(y.range)] <- 10
              } else {
                y.range <- c(0.2, 4)
              }
            }

            plot(dates[date.idxs],
                 y,
                 type = "n",
                 bty  = "l",
                 ylab = "Fitted", 
                 main = paste(j, i),
                 ylim = y.range)

            polygon(c(dates[date.idxs], rev(dates[date.idxs])),
                    c(pred.lcl, rev(pred.ucl)),
                    col = rgb(1, 0, 0, alpha = 0.1), 
                    border = NA)

            lines(dates[date.idxs],
                  pred,
                  type = "l",
                  col = "red",
                  lty = 1, 
                  lwd = 1.5)

            lines(dates[date.idxs], y, lwd = 1.5)

            if (input$set.prop == "Observed/Expected") {
              abline(h = 1, col = "gray", lty = 2)
            }
          }
      }
    } ,

    width = "auto",
    height = "auto"
  )}

  si <- shiny::selectInput
  
  ui <- shiny::fluidPage(
    shiny::titlePanel(paste0('Data through ', last.date.format)), 
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        si("set.prop",
           "Proportion of ED visits or count:",
           choice=c('Proportion','Counts','Observed/Expected'),
           selected ="Proportion" ),

        si("set.borough",
           "Geographic unit:",
           choice=counties.to.test,
           selected ="Citywide" ),

        si("set.syndrome",
           "Syndrome:",
           choice=syndromes,
           selected ="ili"),

        shiny::sliderInput(
          'display.dates',
          'Earliest date to display',
          min=min(dates),
          max=dates[length(dates)-2],
          step=7,
          value=dates[length(dates) - round(length(dates)/5)]
        ),

        shiny::checkboxInput("set.axis",
                             "Uniform axis for all plots?:",
                             value=F),

        si("arrange.plots",
           "Arrange plots by:",
           choice=c('Age','Region'),
           selected ="Age"),

        si("set.ages",
           "Ages:",
           choice=ages.to.test,
           selected =ages.to.test,
           multiple=T)
      ),

      shiny::mainPanel(
        shiny::plotOutput("countyPlot"),
        shiny::column(
          8,
          align = 'justify',
          shiny::hr(),
          shiny::span("The black line shows the observed number of ED visits per day in the indicated stratum, and the red lines denote the mean and 95% prediction intervals for a model adjusting for seasonality, influenza activity, and RSV activity"),
          shiny::hr(),
          shiny::span("This app and package were developed by The Public Health Modeling Unit and The Weinberger Lab at Yale School of Public Health. Contributors include Dan Weinberger, Alyssa Amick, Forrest Crawford, Kelsie Cassell, Marcus Russi, Ernest Asare, Yu-Han Kao. Underlying analysis code can be found at https://github.com/weinbergerlab/ExcessILI")
          )
        )
      )
    )
  shiny::shinyApp(ui, server)
}
