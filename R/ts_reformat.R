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
#' @param covs A character vector. Which, if any, variables in \code{ds} should
#'   be treated as covariates in fitting the baseline model? Default is to not
#'   consider any variables in \code{ds} to be covariates.
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
           covs=character(),
           syndromes,
           resolution='day',
           remove.final=F) {

  att       <- assertthat::assert_that
  is.string <- assertthat::is.string

  att(is.data.frame(line.list))
  att(is.string(datevar))
  att(is.string(statevar))
  att(is.string(sub.statevar))
  att(is.string(agevar))
  att(is.character(covs))

  if (length(covs) > 0)
    att(all(covs %in% names(line.list)))

  att(is.character(syndromes))
  att(is.string(resolution))
  att(any(resolution %in% c('day', 'week', 'month')))
  att(is.logical(remove.final))
  
  ds1 <- line.list
  
  # Parse dates into Date objects, and floor to the nearest day
  ds1[, datevar] <- as.Date(ds1[,datevar])
  ds1[, datevar] <- lubridate::floor_date(ds1[, datevar], unit=resolution)
  
  ds1$all.visits <- 1
  
  if(!(sub.statevar %in% names(ds1))){
    ds1$sub.statevar <- ds1[,statevar]
    sub.statevar <-'sub.statevar'
  }
  
  id_vars       <- c(agevar, datevar, statevar, sub.statevar)
  included_vars <- c(id_vars, syndromes, 'all.visits', covs)

  ds1.molten <- reshape2::melt(ds1[, included_vars], id.vars = id_vars)

  last.date <- max(ds1.molten[,datevar])
  
  # remove last day from the dataset,assuming it is incomplete
  if (remove.final)
    ds1.molten <- filter(ds1.molten, datevar < last.date)
  
  as.formula(
    paste0( paste(id_vars, collapse=" + "), ' ~ ', 'variable' )
  ) -> form1

  ds1.casted <- reshape2::dcast(ds1.molten, form1, fun.aggregate = sum)

  return(ds1.casted)
}
