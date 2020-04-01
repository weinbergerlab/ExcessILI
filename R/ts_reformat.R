#' Format line list data into time series
#'
#' \code{ts_format} takes a line list of case data and formats it into weekly 
#'   or daily time series, which can be used to fit a seasonal baseline.
#'
#' @param line.list A dataframe containing one line for each case (e.g., ED
#'   visit, hospitalization).  At a minimum, each row should have the date of
#'   the visit (\code{YYYY-MM-DD}), the state code (e.g. \code{"NY"}), and a
#'   0/1 variable for each syndrome of interest (e.g. influenza-like illness,
#'   fever, cough). All visits should be included in the dataframe, even if the
#'   case did not have any of the syndromes of interest. For instance, for
#'   emergency department data, every ED visit should have a line represented
#'   in the dataframe
#' 
#' @param datevar A string. What variable contains the date?
#'
#' @param statevar A string. What variable contains the 2-digit state code
#'   (e.g., \code{"NY"})?
#'
#' @param sub.statevar A string. What variable contains the local geography
#'   identifier? (e.g., county, borough)
#'
#' @param agevar A string. What variable contains the age group? Use 'none'
#'   if there is no age grouping in the data
#'
#' @param syndromes A character vector. Which variables contain counts of
#'   syndromic data? (e.g., \code{c('ili', 'respiratory')})
#'
#' @param resolution One of \code{c("day", "week", "month")}. What is the data
#'   binned by?
#'
#' @param remove.final A logical scalar. Remove the final date in the dataset?
#'   This is someties helpful if the data from the last date is unfinalized
#'   or otherwise untrustworthy.
#'
#' @return A dataframe in the "long" format, with a row for each time period
#'   (as in, week or day), and location (e.g. state, county), and age category.
#'   There is a column for date, age category, location, and the number of
#'   counts for each of the selected syndromes. There is also a column that
#'   tallies all visits, regardless of cause
#'
#' @examples
#'  n.obs <- 10000
#'  set.seed(42)
#'
#'  simulated_data <-
#'    as.data.frame(matrix(NA, nrow=n.obs, ncol=5))
#'
#'  names(simulated_data) <- c('state','date','agegrp','ili','resp')
#'
#'  simulated_data$state<- c( rep('CT', times=n.obs*0.3),
#'                  rep("NY", times=n.obs*0.7) )
#'
#'  simulated_data$agegrp <- sample(1:5, n.obs, replace=T)
#'  simulated_data$date   <- 
#'    sample(seq.Date(from=as.Date('2019-01-01'), by='day', length.out=500),                                      
#'           1000,                                                 
#'           replace=T)                                            
#'
#'  simulated_data$ili  <- rbinom(n=n.obs, size=1, prob=0.05)
#'  simulated_data$resp <- rbinom(n=n.obs, size=1, prob=0.1)
#'
#'  ts1 <- ts_format(line.list=simulated_data,
#'                   datevar='date',
#'                   agevar='agegrp',
#'                   statevar='state',
#'                   syndromes=c('ili','resp'))
#'
#' @export
ts_format <-
  function(line.list,
           datevar, statevar, sub.statevar='none', agevar='none',
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
    ds1$sub.statevar <- ds1[,statevar]
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
