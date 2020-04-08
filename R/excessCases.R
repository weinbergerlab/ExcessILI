#' Fits a baseline to case data
#'
#' \code{excessCases} takes a time series of cases (daily or weekly) and fits a
#'   harmonic baseline. There is also an option to import influenza data from
#'   the CDC's NREVSS database and match it by state, or import Google search
#'   queries for RSV for the respective state. Dummy variables adjust for
#'   variations in average incidence between years, and interactions between
#'   RSV or flu allow these effects to vary over time.
#'
#' @param ds A data.frame, with a format similar to the one produced by the
#'   \code{\link{ts_format}} function. There should be a row for each time
#'   period (week or day), location (e.g. state, county), and age category.
#'   There must be a column for date (\code{YYYY-MM-DD}), age category,
#'   location, and the number of counts for each of the selected syndromes.
#'   There should also be a column with a denominator (e.g., total number of ED
#'   visits). If there is no denominator, create a variable with vector of 1s
#'   as a substitute.
#'  
#' @param sub.statevar A string. Which variable in the input data frame
#'   contains the local geography identifier (e.g., county, borough)
#' 
#' @param statevar A string. Which variable in the input data frame contains
#'   the state (2-digit state; e.g. \code{'NY'})?
#'
#' @param agevar A string. Which variable in the input data frame contains the
#'   age group? Use \code{'none'} if there is no age grouping in the data
#'
#' @param datevar A string. Which variable in the input data frame contains
#'   the date?
#'
#' @param covs A character vector. Which, if any, variables in \code{ds} should
#'   be treated as covariates in fitting the baseline model? Default is to not
#'   consider any variables in \code{ds} to be covariates.
#'
#' @param use.syndromes A vector with the variable names for syndromes to be
#'   tested (e.g., \code{c('ILI','respiratory')} ).
#'
#' @param denom.var A string. Which variable on the input dataframe should be
#'   used as the denominator? For instance, all ED visits.
#'
#' @param flu.import A logical scalar. Import the latest influenza testing data
#'   from the CDC NREVSS system? If TRUE, the data will be downloaded and merge
#'   with the input dataframe by state and week. the flu variable will be
#'   included in the regressionwhen fitting the baseline
#'
#' @param rsv.import A logical scalar. Import weekly search volume for 'RSV'
#'   for the states in the input dataframe? This option can only be used if
#'   there are 5 or fewer states on the input dataset. This variable is
#'   included in the regression model when fitting the seasonal baseline.
#'
#' @param adj.flu A logical scalar. How should influenza be adjusted for when
#'   fitting the seasonal baseline? Possible values are \code{'none'} for no
#'   adjustment (default); 'auto': automatically downloads NREVSS data from CDC
#'   and matches by state and week; or specify the name of a variable in the
#'   the input dataframe that contains a variable for influenza. \strong{The
#'   package \code{cdcfluview} is required to use \code{'auto'}.}
#'
#' @param adj.rsv A string. How should RSV be adjusted for when fitting the
#'   seasonal baseline? Possible values are \code{'none'} for no adjustment
#'   (default); \code{'auto'}: automatically downloads the weekly volume of
#'   search queries for 'RSV' for the last 5 years from Google trends and
#'   matches by state. Note that a maximum of 5 states can be included on the
#'   input dataset when using the 'auto option'; Or specify the name of a
#'   variable in the the input dataframe that contains a variable for
#'   influenza. 
#'
#' @param time.res One of \code{c("day", "week", "month")}. What is the data
#'   binned by?
#'
#' @param extrapolation.date The model is fit up to this date, and then
#'   extrapolated for all future dates. Defaults to \code{"2020-03-01"}
#'
#' @return A list of lists with an entry for each syndrome, and sub-lists by
#'   age group and geography:
#'
#' \code{date}: vector of dates used in the model. Use the helper function
#'   \code{\link{excessExtract}} to pull out specific components and organize
#'   them into an array
#'
#' \code{y}: array of observed values for the syndromes
#'
#' \code{resid1}: Observed/model fitted values
#'
#' \code{upi}: upper 95% prediction interval of the fitted value
#'
#' \code{lpi}: lower 95% prediction interval of the fitted value
#'
#' \code{sqrt.rsv}: RSV variale used in the model (if any)
#'
#' \code{log.flu}: flu variable used in the model (if any)
#'
#' \code{unexplained.cases}: observed-expected(fitted)
#'
#' \code{denom}: denominator used in the model
#' 
#' \code{pred.var}: variance of the prediction interval
#'
#' @examples
#'  library(cdcfluview)
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
           covs=character(),
           use.syndromes,
           denom.var,
           adj.flu='none',
           adj.rsv='none',
           time.res='day',
           extrapolation.date='2020-03-01') {

  att       <- assertthat::assert_that
  is.string <- assertthat::is.string

  att(is.data.frame(ds))
  att(is.string(sub.statevar))
  att(is.string(statevar))
  att(is.string(datevar))
  att(is.character(covs))

  if (length(covs) > 0)
    att(all(covs %in% names(ds)))

  att(is.character(use.syndromes))
  att(is.string(denom.var))
  att(is.string(adj.flu))
  att(is.string(adj.rsv))
  att(is.string(time.res) && time.res %in% c('day', 'week', 'month'))
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

  if(adj.flu=='none'){
    ds1.df$flu.var<-1
  }

  if(adj.rsv=='none'){
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
      covs,
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

  combo2.sub <- ds1.df[, cols_of_interest]
  
  if(time.res == 'week'){
    combo2.sub[,datevar] <-
      lubridate::floor_date(combo2.sub[,datevar], unit='week')
  }
  
  ds2 <-
    reshape_ds(ds2          = combo2.sub,
               sub.statevar = sub.statevar,
               agevar       = agevar,
               datevar      = datevar)
 
  ages <- dimnames(ds2)[[3]]
  geos <- dimnames(ds2)[[2]]

  all.glm.results <- pbapply::pblapply(use.syndromes, function(x){
    f <- function(y) {
      model.result <- lapply(
        geos, glm.func,
        ds        = ds2,
        age.test  = y,
        syndrome  = x,
        adj.flu   = adj.flu,
        adj.rsv   = adj.rsv,
        covs      = covs,
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
