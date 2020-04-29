## Import auxillary data Pulls in google searches for RSV in CT can get 5 years
## of historical data at weekly reslution (switches to monthly for >5 years)
rsv.google.import <- function(geo.select) {
  
  rsv.down <-
    gtrendsR::gtrends(keyword = "RSV",
                      geo = paste0("US-", geo.select),
                      time = "today+5-y",
                      gprop = c("web"),
                      category = 0,
                      low_search_volume = FALSE)
  
  rsv <- rsv.down$interest_over_time[, c("date", "hits", "geo")]
  
  rsv$geo  <- gsub("US-", "", rsv$geo, fixed = T)
  rsv$date <- as.Date(rsv$date)
  
  mmwr.week.rsv <- MMWRweek::MMWRweek(rsv$date)[, c("MMWRyear", "MMWRweek")]
  
  rsv <- cbind.data.frame(rsv, mmwr.week.rsv)
  
  names(rsv) <- c("date", "rsv.var", "state", "MMWRyear", "MMWRweek")
  
  return(rsv)
}

# Pull NREVSS testing data (% positive)
nrevss_flu_import <- function() {
  
  if (!requireNamespace("cdcfluview", quietly = TRUE)) {
    stop("Package \"cdcfluview\" needed for this function to work. Please install it.",
      call. = FALSE)
  }

  nrevvs.state <- cdcfluview::who_nrevss(region = c("state"))
  
  clin <- nrevvs.state[["clinical_labs"]]
  
  data(hhs_regions)
  
  cw.file <- hhs_regions
  
  clin2 <- merge(clin, cw.file,
                 by.x = "region",
                 by.y = "state_or_territory")
  
  clin2.subsetvars <- 
    c('region', 'region_number',
      'year', 'week', 'wk_date',
      'total_a','total_b',
      'total_specimens')
  
  clin2 <- clin2[, clin2.subsetvars]
  
  names(clin2)[1:2] <- c("state", "hhs_region")
  
  nrevvs_hhs <- cdcfluview::who_nrevss(region = c("hhs"))
  
  clin.hhs <- nrevvs_hhs[["clinical_labs"]]
  clin.hhs.subsetvars <-
    c('region',
      'wk_date',
      "total_a",'total_b',
      'total_specimens')
  
  clin.hhs <- clin.hhs[, clin.hhs.subsetvars]
  clin.hhs$region <- as.numeric(gsub("Region ", "", clin.hhs$region))
  
  names(clin.hhs) <-
    c("hhs_region",
      "wk_date",
      "hhs_total_a",'hhs_total_b',
      'hhs_total_specimens')
  
  clin3 <- merge(clin2, clin.hhs,
                 by = c("hhs_region", "wk_date"))
  
  clin3$total_a[is.na(clin3$total_a)] <-
    clin3$hhs_total_a[is.na(clin3$total_a)]
  
  clin3$total_b[is.na(clin3$total_b)] <-
    clin3$hhs_total_b[is.na(clin3$total_b)]
  
  clin3$total_specimens[is.na(clin3$total_specimens)] <-
    clin3$hhs_total_specimens[is.na(clin3$total_specimens)]
  
  clin3$state.abb <- state.abb[match(clin3$state, state.name)]
  
  names(clin3) <-
    c("hhs_region",
      "wk_date",
      "state_name",
      "MMWRyear", "MMWRweek",
      "total_a",'total_b',
      'total_specimens',
      'total_a_hhs', "total_b_hhs",
      'total_specimens_hhs',
      "state")
  
  clin3$total_a         <- as.numeric(clin3$total_a)
  clin3$total_b         <- as.numeric(clin3$total_b)
  clin3$total_specimens <- as.numeric(clin3$total_specimens)
  clin3$flu_pct_adj     <- (clin3$total_a + clin3$total_b + 0.5) / 
    (clin3$total_specimens + 0.5)
  clin3$fluN            <- clin3$total_a + clin3$total_b + 0.5
  clin3$flu.var         <- clin3$flu_pct_adj
  
  return(clin3)
}

reshape_ds <- function(ds2, agevar, datevar, sub.statevar) {
  ili.m <- reshape2::melt(
    ds2, 
    id.vars = c(sub.statevar, agevar, datevar, "MMWRyear", "MMWRweek")
  )
  
  casting_formula <-
    as.formula(paste(paste0(datevar, "+MMWRyear+MMWRweek"),
                     sub.statevar, agevar, "variable",
                     sep = "~"))
  
  ili.a <- reshape2::acast(ili.m, casting_formula, fun.aggregate = sum)
  
  dimnames(ili.a)[[1]] <- substr(dimnames(ili.a)[[1]], 1, 10)
  
  return(ili.a)
}


## Evaluate results after controlling for flu and RSV
#' @importFrom magrittr %>%
glm.func <- function(ds, x.test, age.test, denom.var, syndrome, time.res,
                     extrapolation.date,sum.dates, adj.flu, adj.rsv, covs=character(), model.type, seedN,
                     stage1.samples, stage2.samples)
{
  date.string       <- as.Date(dimnames(ds)[[1]])
  month             <- lubridate::month(date.string)
  epiyr             <- lubridate::year(date.string)
  epiyr[month <= 6] <- epiyr[month <= 6] - 1
  epiyr.index       <- epiyr - min(epiyr) + 1
  weekN             <- MMWRweek::MMWRweek(date.string)[, "MMWRweek"]
  day.of.year       <- lubridate::yday(date.string)
  day.of.week       <- as.factor(weekdays(date.string))
  
  clean.array.citywide <- ds[, x.test, , , drop = F]
  epiyr.index.f        <- as.factor(epiyr.index)
  epiyr.index.f2       <- as.factor(epiyr.index)
  y.age                <- t(clean.array.citywide[, , age.test, syndrome])
  n.dates              <- length(y.age)
  y.age.fit            <- y.age[1, ]
  
  # Extrapolate last 1 months
  y.age.fit[date.string >= as.Date(extrapolation.date)] <- NA  
  
  # Same for all ages and boroughs
  sqrt.rsv <- sqrt(clean.array.citywide[, , age.test, "rsv.var"])  

  # Extract the columns of the covariates.
  covs_for_glm <-
    purrr::map(covs, ~clean.array.citywide[, , age.test, .]) %>% setNames(covs)
  
  # Get rid of any zeroes in the data, to play nice with logarithms
  if( min(clean.array.citywide[, , age.test, "flu.var"], na.rm=T) == 0 ){
    clean.array.citywide[, , age.test, "flu.var"] <-
      clean.array.citywide[, , age.test, "flu.var"] + 0.001
  }
  
  # Same for all ages and boroughs
  flu <- clean.array.citywide[, , age.test, "flu.var"]
  
  cont.correct <- 0

  if (min(flu, na.rm = T) == 0)
    cont.correct <- min(flu[flu != 0], na.rm = T)/2 

  log.flu = log(flu + cont.correct)  
  
  # Fill in missing observations for flu at end of TS with most recent
  # observed values
  log.flu <- zoo::na.locf(log.flu, na.rm = F)  
  
  # Fill in missing observations for RSV at end of TS with most recent
  # observed values
  sqrt.rsv <- zoo::na.locf(sqrt.rsv, na.rm = F)  
  
  # Fill in missing observations for covs at end of TS with most recent
  # observed values. NOTE: this may not be appropriate!
  covs_for_glm <- purrr::map(covs_for_glm, zoo::na.locf, na.rm = F)

  #t2 <- 1:length(y.age)
  if(time.res %in% c('day', 'week')){
    t2 <- as.vector( difftime(date.string,min(date.string), units=time.res))
  }else if(time.res=='month'){
    t2 <- round(as.vector(difftime(date.string,min(date.string),units='day'))/30)
  }
  
  if (time.res == "day") {
    period = 365.25
  } else if (time.res == "week") {
    period = 52.1775
  }  else if (time.res =='month'){ 
    period=12 
  }
  
  sin1 <- sin(2*pi * t2/period)
  cos1 <- cos(2*pi * t2/period)
  sin2 <- sin(2*pi * t2 * 2/period)
  cos2 <- cos(2*pi * t2 * 2/period)
  sin3 <- sin(2*pi * t2 * 3/period)
  cos3 <- cos(2*pi * t2 * 3/period)
  
  log.offset <- log(clean.array.citywide[, , age.test, denom.var] + 0.5)
  
  vars_for_glm <-
    list(
      date          = date.string,
      day.of.year   = day.of.year,
      y.age         = y.age[1, ],
      y.age.fit     = y.age.fit,
      sqrt.rsv      = sqrt.rsv,
      log.flu       = log.flu,
      day.of.week   = day.of.week,
      t2            = t2,
      epiyr.index.f = epiyr.index.f, 
      log.offset    = log.offset,
      sin1=sin1, sin2=sin2, cos1=cos1, cos2=cos2, sin3=sin3, cos3=cos3
    )
  
  # Splice together the two lists of varialbes that we would like to have in
  # the dataframe that ultimately passed to 'glm'.
  ds.glm <- do.call(cbind.data.frame, purrr::splice(vars_for_glm, covs_for_glm))

  ds.glm <- ds.glm[!is.na(ds.glm$sqrt.rsv) & !is.na(ds.glm$log.flu), ]
  ds.glm$epiyr.index.f <- factor(ds.glm$epiyr.index.f)

  ###########################################
  # Model specification
  ###########################################

  # Term 1a,1b adjust for flu and RSV. If both of them aren't specified, 
  # then only one 'epiyr.index.f' is to be specified.
  covars.term1a <- ifelse(adj.flu == 'none', 'epiyr.index.f', 'epiyr.index.f*log.flu')
  covars.term1b <- ifelse(adj.rsv == 'none', 'epiyr.index.f', 'epiyr.index.f*sqrt.rsv')

  covars.term1  <- ifelse(identical(covars.term1a, covars.term1b),
                          covars.term1a, 
                          paste(covars.term1a, covars.term1b, sep=" + "))

  # If the data has day-resolution, use it as one of the mock variables
  covars.term2  <- ifelse(time.res == 'day', 'day.of.week', NA)

  # Sinusoidals
  covars.term3  <- paste("sin1", "cos1", "sin2", "cos2","sin3","cos3" , sep=" + ")

  # Add in user-specified covariates, if available
  covars.term4  <- ifelse(length(covs > 0), paste(covs, collapse = " + "), NA)

  offset.term <- "offset(log.offset)"
  
  # Concatenate everything together as a string
  covars <- c(covars.term1, covars.term2, covars.term3, covars.term4, offset.term)

  covars_str <- paste(covars[!is.na(covars)], collapse=" + ")
  
  # Rsv effect varies by epiyr
  form1 <- as.formula(paste0("y.age.fit ~ ", covars_str))
  form2 <- as.formula(paste0("y.age ~ ",     covars_str))
 
  ###########################################
  # Model fitting
  ###########################################

  if (sum(ds.glm$y.age, na.rm=T) >= 100) {
    
    if(model.type=='poisson'){
    mod1 <- glm(form1,
                data = ds.glm,
                family = poisson(link = "log"))
    }else{
    mod1 <- MASS::glm.nb(form1,
                     data = ds.glm)
    }
    # 500 samples total
    coef1               <- coef(mod1)
    coef1[is.na(coef1)] <- 0
    
    v.cov.mat                   <- vcov(mod1)
    v.cov.mat[is.na(v.cov.mat)] <- 0
    
    disp <- mod1$deviance/mod1$df.residual
    
    set.seed (seedN)
    pred.coefs.reg.mean <-
      MASS::mvrnorm(n = stage1.samples,
                    mu = coef1,
                    Sigma = v.cov.mat)
    
    mod.mat.pred <- model.matrix(form2, data = ds.glm, family = "poisson")
    
    preds.stage1.regmean <- mod.mat.pred %*% t(pred.coefs.reg.mean)
    preds.stage1.regmean <- apply(preds.stage1.regmean, 2,
                                  function(x) x + ds.glm$log.offset)
    
    if(model.type=='poisson'){
    preds.stage2 <- rpois(n = length(preds.stage1.regmean) * stage2.samples,
                          exp(preds.stage1.regmean))
    }else{
    preds.stage2 <- rnbinom(n = length(preds.stage1.regmean) * stage2.samples,
                            size = mod1$theta, mu = exp(preds.stage1.regmean))
    }
    preds.stage2 <- matrix(preds.stage2,
                           nrow = nrow(preds.stage1.regmean),
                           ncol = ncol(preds.stage1.regmean) * stage2.samples)
    
    preds.stage2.q <-
      t(apply(preds.stage2, 1,
              quantile,
              probs = c(0.025, 0.5, 0.975)))
    
    preds.stage2.var <-
      apply(preds.stage2, 1,
              var, na.rm=T)
    
    eval.indices <- 
      which(ds.glm$date >= sum.dates)
    
    sum.pred.iter <- apply(preds.stage2[eval.indices,], 2,
          sum, na.rm=T)

    sum.obs <- ds.glm$y.age[eval.indices]
    
    resid1 <- log( (ds.glm$y.age + 0.5) / 
                     (preds.stage2.q[, "50%"] + 0.5))
    
    unexplained.cases <- ds.glm$y.age - preds.stage2.q[, "50%"]
    
    out.list <-
      list(date              = ds.glm$date,
           y                 = ds.glm$y.age,
           pred              = preds.stage2.q[, "50%"],
           resid1            = resid1,
           upi               = preds.stage2.q[, "97.5%"],
           lpi               = preds.stage2.q[, "2.5%"],
           sqrt.rsv          = ds.glm$sqrt.rsv,
           log.flu           = ds.glm$log.flu,
           unexplained.cases = unexplained.cases, 
           denom             = exp(ds.glm$log.offset),
           pred.var          = preds.stage2.var,
           sum.pred.iter     = sum.pred.iter,
           pred.iter         = preds.stage2,
           sum.obs           = sum.obs,
           disp               = disp,
           sparse.group      = F)
  } else {
    out.list <-
      list(date              = ds.glm$date,
           y                 = ds.glm$y.age,
           pred              = NA,
           resid1            = NA,
           upi               = NA,
           lpi               = NA,
           sqrt.rsv          = ds.glm$sqrt.rsv,
           log.flu           = ds.glm$log.flu,
           unexplained.cases = NA, 
           denom             = exp(ds.glm$log.offset),
           pred.var          = NA,
           sum.pred.iter     = NA,
           pred.iter         = NA,
           sum.obs           = sum.obs,
           disp              = NA,
           sparse.grp        = T)
  }  
  return(out.list)
}
