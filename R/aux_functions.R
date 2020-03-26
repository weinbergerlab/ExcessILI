## Import auxillary data Pulls in google searches for RSV in CT can get 5 years
## of historical data at weekly reslution (switches to monthly for >5 years)
rsv.google.import <- function(geo.select) {

    rsv.down <- gtrends(keyword = "RSV",
                        geo = paste0("US-", geo.select),
                        time = "today+5-y",
                        gprop = c("web"),
                        category = 0,
                        low_search_volume = FALSE)

    rsv <- rsv.down$interest_over_time[, c("date", "hits", "geo")]
    rsv$geo <- gsub("US-", "", rsv$geo, fixed = T)
    rsv$date <- as.Date(rsv$date)

    mmwr.week.rsv <- MMWRweek(rsv$date)[, c("MMWRyear", "MMWRweek")]

    rsv <- cbind.data.frame(rsv, mmwr.week.rsv)

    names(rsv) <- c("date", "rsv.var", "state", "MMWRyear", "MMWRweek")

    return(rsv)
}

# Pull NREVSS testing data (% positive)
nrevss_flu_import <- function() {

    nrevvs.state <- who_nrevss(region = c("state"))

    clin <- nrevvs.state[["clinical_labs"]]

    data(hhs_regions)
    cw.file <- hhs_regions

    clin2 <- merge(clin, cw.file, by.x = "region", by.y = "state_or_territory")
    clin2 <- clin2[, c("region", "region_number", "year", "week", "wk_date", "percent_positive")]

    names(clin2)[1:2] <- c("state", "hhs_region")

    nrevvs_hhs <- who_nrevss(region = c("hhs"))
    clin.hhs <- nrevvs_hhs[["clinical_labs"]]
    clin.hhs <- clin.hhs[, c("region", "wk_date", "percent_positive")]
    clin.hhs$region <- as.numeric(gsub("Region ", "", clin.hhs$region))
    names(clin.hhs) <- c("hhs_region", "wk_date", "hhs_percent_positive")

    clin3 <- merge(clin2, clin.hhs, by = c("hhs_region", "wk_date"))

    clin3$percent_positive[is.na(clin3$percent_positive)] <-
      clin3$hhs_percent_positive[is.na(clin3$percent_positive)]

    clin3$state.abb <- state.abb[match(clin3$state, state.name)]
    names(clin3) <- c("hhs_region", "wk_date", "state_name", "MMWRyear",
                      "MMWRweek", "flu_pct_pos", "hhs_percent_positive",
                      "state")
    clin3$flu.var <- as.numeric(as.character(clin3$flu_pct_pos))

    return(clin3)
}

reshape_ds <- function(ds2, agevar, datevar, sub.statevar) {
    ili.m <- melt(ds2, id.vars = c(sub.statevar, agevar, datevar, "MMWRyear", "MMWRweek"))

    casting_formula <-
      as.formula(paste(paste0(datevar, "+MMWRyear+MMWRweek"),
                       sub.statevar, agevar, "variable",
                       sep = "~"))

    ili.a <- acast(ili.m, casting_formula, fun.aggregate = sum)
    dimnames(ili.a)[[1]] <- substr(dimnames(ili.a)[[1]], 1, 10)

    return(ili.a)
}


## Evaluate results after controlling for flu and RSV
glm.func <- function(ds, x.test, age.test, denom.var, syndrome, time.res,extrapolation.date) {
    date.string <- as.Date(dimnames(ds)[[1]])
    month <- month(date.string)
    epiyr <- year(date.string)
    epiyr[month <= 6] <- epiyr[month <= 6] - 1
    epiyr.index <- epiyr - min(epiyr) + 1
    weekN <- MMWRweek(date.string)[, "MMWRweek"]
    day.of.year <- yday(date.string)
    day.of.week <- as.factor(weekdays(date.string))
    
    clean.array.citywide <- ds[, x.test, , , drop = F]
    epiyr.index.f <- as.factor(epiyr.index)
    epiyr.index.f2 <- as.factor(epiyr.index)
    y.age = t(clean.array.citywide[, , age.test, syndrome])
    n.dates <- length(y.age)
    y.age.fit <- y.age[1, ]
    y.age.fit[date.string >= as.Date("2020-03-01")] <- NA  #extrapolate last 1 months
    sqrt.rsv = sqrt(clean.array.citywide[, , age.test, "rsv.var"])  #same for all ages and boroughs
    sqrt.flu = sqrt(clean.array.citywide[, , age.test, "flu.var"])  #same for all ages and boroughs
    sqrt.flu <- na.locf(sqrt.flu, na.rm = F)  #fill in missing observations for flu at end of TS with most recent observed values
    sqrt.rsv <- na.locf(sqrt.rsv, na.rm = F)  #fill in missing observations for RSV at end of TS with most recent observed values
    t2 <- 1:length(y.age)
    if (time.res == "day") {
        period = 365.25
    } else if (time.res == "week") {
        period = 52.1775
    }
    sin1 <- sin(2 * pi * t2/period)
    cos1 <- cos(2 * pi * t2/period)
    sin2 <- sin(2 * pi * t2 * 2/period)
    cos2 <- cos(2 * pi * t2 * 2/period)
    sin3 <- sin(2 * pi * t2 * 3/period)
    cos3 <- cos(2 * pi * t2 * 3/period)
    
    log.offset <- log(clean.array.citywide[, , age.test, denom.var])
    ds.glm <-
      cbind.data.frame(
        date = date.string,
        day.of.year = day.of.year,
        y.age = y.age[1, ],
        y.age.fit = y.age.fit,
        sqrt.rsv,
        sqrt.flu,
        day.of.week,
        t2,
        epiyr.index.f, 
        sin1, sin2, cos1, cos2, sin3, cos3,
        log.offset)

    ds.glm <- ds.glm[!is.na(ds.glm$sqrt.rsv) & !is.na(ds.glm$sqrt.flu), ]
    ds.glm$epiyr.index.f <- factor(ds.glm$epiyr.index.f)
    # ds.glm<-ds.glm[complete.cases(ds.glm),]

    if (time.res == "day") {
        # rsv effect varies by epiyr
        form1 <- as.formula(paste0("y.age.fit ~",
                                   "epiyr.index.f*sqrt.rsv +",
                                   "epiyr.index.f*sqrt.flu +", # flu effect, varies by epiyear
                                   "day.of.week +",
                                   "sin1+cos1 + sin2+cos2 + sin3+cos3"))
    
        # rsv effect varies by epiyr
        form2 <- as.formula(paste0("y.age ~",
                                   "epiyr.index.f*sqrt.rsv + ",
                                   "epiyr.index.f*sqrt.flu + ", # flu effect, varies by epiyear
                                   "day.of.week + ",
                                   "sin1+cos1 + sin2+cos2 + sin3+cos3"))
    } else {
        #rsv effect varies by epiyr
        form1 <- as.formula(paste0("y.age.fit ~",
                                   "epiyr.index.f*sqrt.rsv + ",
                                   "epiyr.index.f*sqrt.flu + ", # flu effect, varies by epiyear
                                   "sin1+cos1 + sin2+cos2 + sin3+cos3 "))

  #rsv effect varies by epiyr
        form2 <- as.formula(paste0("y.age ~",
                                   "epiyr.index.f*sqrt.rsv + ",
                                   "epiyr.index.f*sqrt.flu + ", #flu effect, varies by epiyear
                                   "sin1+cos1 + sin2+cos2 + sin3+cos3"))
    }

    mod1 <- glm(form1, data = ds.glm, family = poisson(link = "log"), offset = log.offset)

    # 500 samples total
    coef1 <- coef(mod1)
    coef1[is.na(coef1)] <- 0
    v.cov.mat <- vcov(mod1)
    v.cov.mat[is.na(v.cov.mat)] <- 0
    pred.coefs.reg.mean <- mvrnorm(n = 100, mu = coef1, Sigma = v.cov.mat)
    mod.mat.pred <- model.matrix(form2, data = ds.glm, family = "poisson")
    preds.stage1.regmean <- mod.mat.pred %*% t(pred.coefs.reg.mean)
    preds.stage1.regmean <- apply(preds.stage1.regmean, 2,
                                  function(x) x + ds.glm$log.offset)
    preds.stage2 <- rpois(n = length(preds.stage1.regmean) * 5,
                          exp(preds.stage1.regmean))
    preds.stage2 <- matrix(preds.stage2,
                           nrow = nrow(preds.stage1.regmean),
                           ncol = ncol(preds.stage1.regmean) * 5)
    preds.stage2.q <- t(apply(preds.stage2, 1, quantile, probs = c(0.025, 0.5, 0.975)))
    
    resid1 <- log((ds.glm$y.age + 0.5)/(preds.stage2.q[, "50%"] + 0.5))
    unexplained.cases <- ds.glm$y.age - preds.stage2.q[, "50%"]
    
    out.list <-
      list(date = ds.glm$date,
           y = ds.glm$y.age,
           pred = preds.stage2.q[, "50%"],
           resid1 = resid1,
           upi = preds.stage2.q[, "97.5%"],
           lpi = preds.stage2.q[, "2.5%"],
           sqrt.rsv = ds.glm$sqrt.rsv,
           sqrt.flu = ds.glm$sqrt.flu,
           unexplained.cases = unexplained.cases, 
           denom = exp(ds.glm$log.offset))

    return(out.list)
}

