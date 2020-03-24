#Format line list data into time series
ts_format<-function(ds1, datevar,geovar, agevar, syndromes,resolution='day',remove.final=T){
  ds1[, datevar]<-as.Date(ds1$datevar)
  ds1[, datevar]<-floor_date(ds1[, datevar], unit=resolution)
  ds1$all.visits<-1
  ds1.m<-melt(ds2[,c(datevar, geovar,agevar,syndromes,'all.visits', )], id.vars=c(datevar, geovar, agevar))
  #last.date<- max(ds1.m$adate)
  if(remove.final){
  ds1.m<-ds1.m[ds1.m$adate < last.date,] #remove last day from the dataset,assuming it is incomplete
  }
  form1<-as.formula(paste(agevar,datevar, geovar,'variable', sep='~'))
  ds1.c<-acast(ds1.m, form1, fun.aggregate = sum)
  return(ds1.c)
}

excess_cases<-function(ds2){
  syndromes<- dimnames(ds2)[[4]]
  ages <-   dimnames(ds2)[[3]]
  geos<-dimnames(ds2)[[2]]
  all.glm.res<- pblapply(syndromes, function(x){
    ww<- lapply(ages, function(y){
      lapply(geos, glm.func, age.test=y, syndrome=x)
    }
    ) 
    names(ww)<- ages
    return(ww)
  }
  )
  names(all.glm.res)<-syndromes
  return(all.glm.res)
}

syndromic_dashboard<-function(ds){
  counties.to.test<-c("Bronx","Brooklyn", "Manhattan","Queens","Staten Island", "Citywide" )
  syndromes<-names(ds)
  dates<-as.Date(names(ds[[1]][[1]][[1]]$y))
  n.times<-length(dates)
  last.date.format<-max(dates)
  last.date.format<-format(last.date.format,
                           "%b %d, %Y")
  shinyApp(ui, server)
}

