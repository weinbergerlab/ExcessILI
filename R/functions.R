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

excessCases<-function(ds,geovar,statevar='state',agevar, datevar, use.syndromes,denom.var, flu.import=T, rsv.import=T, adj.flu=T, adj.rsv=T){
  ds<-as.data.frame(ds)
  ds[,datevar]<-as.Date(ds[,datevar])
  mmwr.date<-MMWRweek(ds[,datevar])
  ds1.df<-cbind.data.frame(ds,mmwr.date)
  
  if(is.null(geovar)){
    ds1.df$geovar<-ds1.df[,statevar]
  }
  if(rsv.import){
    #TODO:API only allows pulling 5 or fewer states at a time
  rsv<-rsv.google.import(geo.select=unique(ds1.df[,statevar]))
  ds1.df<-merge(ds1.df, rsv, by.x=c('MMWRyear','MMWRweek', statevar),by.y=c('MMWRyear','MMWRweek', 'state'), all.x=T)
  }
  if(flu.import){
    flu<-nrevss_flu_import()
    ds1.df<-merge(ds1.df, flu, by=c('MMWRyear','MMWRweek',statevar), all.x=T)
  }
  if(adj.flu==F){
    ds1.df$flu.var<-1
  }
  if(adj.rsv==F){
    ds1.df$rsv.var<-1
  }
  if(is.null(denom.var)){
    ds1.df$denom <-1
    denom.var <-'denom'
  }
  combo2.sub<-ds1.df[, c(agevar, datevar,'MMWRyear', 'MMWRweek', geovar, use.syndromes,denom.var, 'flu.var','rsv.var' )]
  ds2<-reshape_ds(ds2=combo2.sub, geovar=geovar, agevar=agevar, datevar=datevar)
 
  ages <-   dimnames(ds2)[[3]]
  geos<-dimnames(ds2)[[2]]
  all.glm.res<- pblapply(use.syndromes, function(x){
    ww<- lapply(ages, function(y){
      lapply(geos, glm.func, ds=ds2,age.test=y, syndrome=x)
    }
    ) 
    names(ww)<- ages
    return(ww)
  }
  )
  names(all.glm.res)<-use.syndromes
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

