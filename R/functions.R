#' Format line list data into time series
#'
#' \code{ts_format} Formats line list data of cases into time series  
#' This function takes a line list of case data and formats it into weekly 
#' or daily time series, which can be used to fit a seasonal baseline
#'
#' @param line.list A dataframe containing one line for each case (e.g., ED visit, hospitalization). 
#' At a minimum, each row should have the date of the visit (YYYY-MM-DD), the state code (e.g. "NY"), and a 0/1 variable for each syndrome of 
#' interest (e.g. influenza-like illness, fever, cough). All visits
#' should be included in the dataframe, even if the case did not have any of the syndromes of interest. For instance,
#' for emergency department data, every ED visit should have a line represented in the dataframe
#' 
#' @param datevar A string. What variable contains the date?
#'
#' @param statevar A string. What variable contains the 2-digit state code (e.g., "NY")?
#'
#' @param sub.statevar A string. What variable contains the local geography identifier (e.g., county, borough)
#'
#' @param agevar A string. What variable contains the age group? Use 'none'
#'   if there is no age grouping in the data
#'
#' @param syndromes A character vector. Which variables contain counts of
#'   syndromic data? (e.g., c('ili', 'respiratory'))
#'
#' @param resolution One of \code{c("day", "week", "month")}. What is the data
#'   binned by?
#'
#' @param remove.final A logical scalar. Remove the final date in the dataset? This is someties 
#' helpful if the data from the last date.
#'
#' @return A dataframe in the "long" format, with a row for each date (week or day), and location (e.g. state, county),
#' and age category. There is a column for date, age category, location, and the number of counts for each of 
#' the selected syndromes. There is also a column that tallies all visits regardless of cause
#' 
#'
#' @examples
#'  n.obs <- 10000
#'  set.seed(42)
#'
#'  simulated_data <-
#'    as.data.frame(matrix(NA, nrow=n.obs, ncol=5))
#'
#'  names(sim1) <- c('state','date','agegrp','ili','resp')
#'
#'  sim1$state<- c( rep('CT', times=n.obs*0.3),
#'                  rep("NY", times=n.obs*0.7) )
#'
#'  sim1$agegrp <- sample(1:5, n.obs, replace=T)
#'  sim1$date   <- sample(seq.Date(from=as.Date('2019-01-01'),
#'                                 by='day',
#'                                 length.out=500), 
#'                        1000,
#'                        replace=T)
#'
#'  sim1$ili  <- rbinom(n=n.obs, size=1, prob=0.05)
#'  sim1$resp <- rbinom(n=n.obs, size=1, prob=0.1)
#'
#'  ts1 <- ts_format(line.list=sim1,
#'                   datevar='date',
#'                   agevar='agegrp',
#'                   statevar='state',
#'                   syndromes=c('ili','resp'))
#'
#' @export
ts_format <-
  function(line.list,
           datevar, statevar, sub.statevar, agevar,
           syndromes,
           resolution='day',
           remove.final=F) {

  is.string <- assertthat::is.string

  assertthat::assert_that(is.data.frame(line.list))
  assertthat::assert_that(is.string(datevar))
  assertthat::assert_that(is.string(statevar))
  assertthat::assert_that(is.string(sub.statevar))
  assertthat::assert_that(is.string(agevar))
  assertthat::assert_that(is.character(syndromes))
  assertthat::assert_that(is.string(resolution))
  assertthat::assert_that(any(resolution %in% c('day', 'week', 'month')))
  assertthat::assert_that(is.logical(remove.final))
  
  ds1 <- line.list
  
  # Parse dates into Date objects, and floor to the nearest day
  ds1[, datevar] <- as.Date(ds1[,datevar])
  ds1[, datevar] <- lubridate::floor_date(ds1[, datevar], unit=resolution)
  
  ds1$all.visits <-1
  
  if(!('sub.statevar' %in% names(ds1))){
    ds1$sub.statevar <- statevar
    sub.statevar <-'sub.statevar'
  }
  
  ds1.m <- reshape2::melt(
    ds1[,   c(datevar, statevar, sub.statevar, agevar,syndromes,'all.visits')],
    id.vars=c(datevar, statevar, sub.statevar, agevar)
  )

  last.date <- max(ds1.m[,datevar])
  
  if(remove.final){
    # remove last day from the dataset,assuming it is incomplete
    ds1.m <- ds1.m[ ds1.m[,datevar] < last.date,] 
  }
  
  form1 <-
    as.formula(
      paste0(
        paste(agevar, datevar, statevar, sub.statevar, sep='+'),
        '~',
        'variable'
      )
    )

  ds1.c <- reshape2::dcast(ds1.m, form1, fun.aggregate = sum)

  return(ds1.c)
}

#' Fits a baseline to the data
#'
#' \code{excessCases} This function takes a time series of cases (daily or weekly) and fits a harmonic baseline
#' There is also an option to import influenza data from the CDC's NREVSS database and match it by state,
#' or import Google search queries for RSV for the respective state. Dummy variables adjust for variations in average
#' incidence between years, and interactions between RSV or flu allow these effects to vary over time.
#'
#' EXTENDED DESCRIPTION, THE FOLLOWING IS AN EXAMPLE: This is a generic
#' function: methods can be defined for it directly or via the
#' \code{\link{Summary}} group generic. For this to work properly, the
#' arguments \code{...} should be unnamed, and dispatch is on the first
#' argument.
#'
#' @param ds A dataframe, with a format similar to the one produced by the ts_format() function. There should be a row for each date (week or day), and location (e.g. state, county),
#' and age category. There must be a column for date (YYYY-MM-DD), age category, location, and the number of counts for each of 
#' the selected syndromes. There should also be a column with a denominator (e.g., total number of ED visits). If there is no denominator, create a vector of 1s.
#'  
#' @param sub.statevar A string. Which variable in the input data frame contains the local geography identifier (e.g., county, borough)
#' 
#' @param statevar A string. Which variable in the input data frame contains the state (2-digit state; e.g. 'NY')?
#'
#' @param agevar A string. Which variable in the input data frame contains the age group? Use 'none'
#'   if there is no age grouping in the data
#'
#' @param datevar A string. Which variables in the input data frame contains the date?
#'
#' @param use.syndromes A vector with the variable names for syndromes to be tested (e.g., c('ILI','respiratory') ).

#' @param denom.var Which variable on the input dataframe should be used as the denominator. For instance, all ED visits. If there is no denominator, use a vector of 1s

#' @param flu.import A logical scalar. Import the latest influenza testing data from the CDC NREVSS system? If TRUE, the data will be downloaded and merge with the input 
#' dataframe by state and week. the flu variable will be included in the regressionwhen fitting the baseline

#' @param rsv.import A logical scalar. Import weekly search volume for 'RSV' for the states in the input dataframe? This option can only be used if there are 5 or fewer states on the input dataset. This variable is included in the regression model when fitting the seasonal baseline.

#' @param adj.flu How should influenza be adjusted for when fittig the seasonal baseline? Possible values are 'none' for no adjustment (default); 'auto': automatically downloads NREVSS data from CDC and matches by state and week; 
#' or specify the name of a variable in the the input dataframe that contains a variable for influenza. Adjust for RSV when fitting the seasonal baseline? This is automatically set to TRUE when

#' @param adj.rsv How should RSv be adjusted for when fitting the seasonal baseline? A string.  Possible values are 'none' for no adjustment (default); 'auto': automatically downloads 
#' the weekly volume of search queries for 'RSV' for the last 5 years from Google trends and matches by state. Note that a maximum of 5 states can be included on the input dataset 
#' when using the 'auto option'; Or specify the name of a variable in the the input dataframe that contains a variable for influenza. 

#' @param time.res One of \code{c("day", "week", "month")}. What is the data
#'   binned by?

#' @param extrapolation.date The model is fit up to this date, and then extrapolated for all future dates.
#'
#' @return A list of lists with an entry for each syndrome, and sub-lists by age group and geography 
#' date: vector of dates used in the model. Use the helper function excessExtract to pull out specific components and organize them into an array
#' y: array of observed values for the syndromes
#' resid1 :Observed/model fitted values
#' upi: upper 95% prediction interval of the fitted value
#' lpi: lower 95% prediction interval of the fitted value
#' sqrt.rsv: RSV variale used in the model (if any)
#' log.flu: flu variable used in the model (if any)
#' unexplained.cases: observed-expected(fitted)
#' denom: denominator used in the model
#'
#' @examples
#'  ili.data <- ilinet(region = c("state"))
#'  ili.data$state <- state.abb[match(ili.data$region, state.name)]
#'
#'  ili.data <- ili.data[, c("state", "week_start", "ilitotal", "total_patients")]
#'  ili.data <- ili.data[ili.data$state %in% c("CA", "NY", "WA", "NJ", "CT"), ]
#'
#'  excess_cases <-
#'    excessCases(ds = ili.data,
#'                datevar = "week_start",
#'                agevar = "none",
#'                statevar = "state",
#'                denom.var = "total_patients",
#'                use.syndromes = c("ilitotal"),
#'                time.res = "week")
#' @export
excessCases <-
  function(ds,
           sub.statevar='none',
           statevar='state',
           agevar='none',
           datevar,
           use.syndromes,
           denom.var,
           adj.flu='none',
           adj.rsv='none',
           time.res='day',
           extrapolation.date='2020-03-01') {

  is.string <- assertthat::is.string

  assertthat::assert_that(is.data.frame(ds))
  assertthat::assert_that(is.string(sub.statevar))
  assertthat::assert_that(is.string(statevar))
  assertthat::assert_that(is.string(datevar))
  assertthat::assert_that(is.character(use.syndromes))
  assertthat::assert_that(is.string(denom.var))
  assertthat::assert_that(is.string(adj.flu))
  assertthat::assert_that(is.string(adj.rsv))
  assertthat::assert_that(is.string(time.res) &&
                          time.res %in% c('day', 'week', 'month'))
  # Need a better test here
  # assertthat::assert_that(!is.null(extrapolation.date)) 

  flu.var='flu.var'
  flu.import <- FALSE
  
  if( adj.flu == 'auto'){
    flu.import <- TRUE
  }
  rsv.var <- 'rsv.var'
  rsv.import <- FALSE
  
  if( adj.rsv == 'auto'){
    rsv.import <- TRUE
  }

  
  if( length(unique(ds[,statevar])) > 5 && identical(rsv.import, T))
    stop('Maximum of 5 states can be used when rsv.import=T')
  
  ds<-as.data.frame(ds)
  ds[,datevar]<-as.Date(ds[,datevar])
  mmwr.date<-MMWRweek::MMWRweek(ds[,datevar])
  ds1.df<-cbind.data.frame(ds,mmwr.date)
  
  if(sub.statevar=='none'){
    sub.statevar<-'sub.statevar'
    ds1.df$sub.statevar<-ds1.df[,statevar]
  }

  if( agevar =='none'){
    ds1.df[,agevar]<-'1'
  }

  if(rsv.import){
    geos <- unique(ds1.df[,statevar])

    if (length(geos) > 5)
      stop(paste0("Cannot query more than 5 geographic regions/states ",
                  "at once using the GTrends api. Your regions:\n",
                  paste(geos, collapse=', ')))

    rsv <-rsv.google.import(geo.select=geos)
    ds1.df <-merge(ds1.df, rsv,
                   by.x=c('MMWRyear','MMWRweek', statevar),
                   by.y=c('MMWRyear','MMWRweek', 'state'),
                   all.x=T)
  }

  if(flu.import){
    flu<-nrevss_flu_import()
    ds1.df<-merge(ds1.df, flu,
                  by=c('MMWRyear','MMWRweek',statevar),
                  all.x=T)
  }

  if(adj.flu=='auto'){
    ds1.df$flu.var<-1
  }

  if(adj.rsv=='auto'){
    ds1.df$rsv.var<-1
  }
  
  if( !(adj.rsv %in% c('auto','none')) ){
    ds1.df$rsv.var <- ds1.df[,adj.rsv]
  } 
  
  if( !(adj.flu %in% c('auto','none')) ){
    ds1.df$flu.var <- ds1.df[,adj.flu]
  } 
  
  if(!exists("denom.var")){
    ds1.df$denom <-1
    denom.var <-'denom'
  }

  cols_of_interest <-
    c(agevar, datevar,
      'MMWRyear', 'MMWRweek',
      sub.statevar,
      use.syndromes,
      denom.var,
      'flu.var',
      'rsv.var')

  if (any( !(cols_of_interest %in% names(ds1.df)) ))
    stop(paste0("Some of 'cols_of_interest' were not in 'ds1.df'.\n\n",
                "'cols_of_interest': ", paste(cols_of_interest, collapse=', '),
                "\n\nnames(ds1.df): ", paste(names(ds1.df), collapse=', '),
                "\n\nmissing: ",
                paste(setdiff(
                        cols_of_interest,
                        intersect(names(ds1.df), cols_of_interest)),
                      collapse=", ")))

  combo2.sub <-
    ds1.df[, c(agevar, datevar,
               'MMWRyear', 'MMWRweek',
               sub.statevar,
               use.syndromes,
               denom.var,
               'flu.var',
               'rsv.var')]
  
  if(time.res=='week'){
    combo2.sub[,datevar] <-
      lubridate::floor_date(combo2.sub[,datevar], unit='week')
  }
  
  ds2 <-
    reshape_ds(ds2=combo2.sub,
               sub.statevar=sub.statevar,
               agevar=agevar,
               datevar=datevar)
 
  ages <- dimnames(ds2)[[3]]
  geos <- dimnames(ds2)[[2]]

  all.glm.results <- pbapply::pblapply(use.syndromes, function(x){
    f <- function(y) {
      model.result <- lapply(
        geos, glm.func,
        ds        = ds2,
        age.test  = y,
        syndrome  = x,
        denom.var = denom.var,
        time.res  = time.res,
        extrapolation.date = extrapolation.date
      )

      names(model.result) <- geos
      return(model.result)
    }

    ww <- lapply(ages, f)
    names(ww) <- ages        

    return(ww)
  })

  names(all.glm.results) <- use.syndromes
  return(all.glm.results)
}

#' Create interactive Shiny app to explore results
#'
#' \code{dashboardPlot} Creates an interactive Shiny plot
#' to explore results generated in the function excessCases. 
#' Drop down menus allow for viewing different syndromes, age groups,
#' and geographies, and for looking at plots of raw counts,
#' proportions, or Observed/Expected
#'

#' @param all.glm.res Provide the object created in the function excessCases
#'
#' @return Launches an interactive shiny app
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
      plucker <- function(var) function(x) sapply(x, "[[", var, simplify='array')
      
      data_to_pluck <- ds[[input$set.syndrome]]

      ili2.resid    <- sapply(data_to_pluck, plucker("resid1"))
      ili2.pred     <- sapply(data_to_pluck, plucker("pred"))
      ili2.pred.lcl <- sapply(data_to_pluck, plucker("lpi"))
      ili2.pred.ucl <- sapply(data_to_pluck, plucker("upi"))
      obs.ili       <- sapply(data_to_pluck, plucker("y"))
      denom         <- sapply(data_to_pluck, plucker("denom"))

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

#' ONE LINER DESCRIPTION
#'
#' \code{excessExtract} Helper function to format the results
#' from excessCases
#'
#' EXTENDED DESCRIPTION This function extracts specific estimates
#' produced in excessCases and aves them in a multidimensional
#' array.
#'
#' @param ds An onject created by the function excessCases
#' @param syndrome Character. For which syndrome should results be extracted? The name for this variable should match a syndrom provided to excessCases
#' @param extract.quantity Which element of the output (e.g., 'lpi') from the output of excessCases do you want to extract. See the values in the help for excessCases to see the options.
#'
#' @return a multideimensional array with dimensions for time, age group, and geography (state, and sub-state)
#'
#' @export
excessExtract <- function(ds, syndrome, extract.quantity) {
    out.ds <-
      sapply(ds[[syndrome]],
             function(x) sapply(x, "[[", extract.quantity), simplify = "array")

    return(out.ds)
}
