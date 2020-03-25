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

excessCases<-function(ds,sub.statevar='none',statevar='state',agevar, datevar, use.syndromes,denom.var, flu.import=T, rsv.import=T, adj.flu=T, adj.rsv=T, flu.var='flu.var', rsv.var='rsv.var', time.res='day'){
  if( length(unique(ds[,statevar]))>5 & rsv.import==T) stop('Maximum of 5 states can be used when rsv.import=T')
  ds<-as.data.frame(ds)
  ds[,datevar]<-as.Date(ds[,datevar])
  mmwr.date<-MMWRweek(ds[,datevar])
  ds1.df<-cbind.data.frame(ds,mmwr.date)
  
  if(sub.statevar=='none'){
    sub.statevar<-'sub.statevar'
    ds1.df$sub.statevar<-ds1.df[,statevar]
  }
  if( !(agevar %in% names(ds1.df))){
    ds1.df[,agevar]<-'1'
  }
  if(rsv.import){
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
  if(!exists("denom.var")){
    ds1.df$denom <-1
    denom.var <-'denom'
  }
  combo2.sub<-ds1.df[, c(agevar, datevar,'MMWRyear', 'MMWRweek', sub.statevar, use.syndromes,denom.var, 'flu.var','rsv.var' )]
  ds2<-reshape_ds(ds2=combo2.sub, sub.statevar=sub.statevar, agevar=agevar, datevar=datevar)
 
  ages <-   dimnames(ds2)[[3]]
  geos<-dimnames(ds2)[[2]]
  all.glm.res<- pblapply(use.syndromes, function(x){
    ww<- lapply(ages, function(y){
      q<-lapply(geos, glm.func, ds=ds2,age.test=y, syndrome=x, denom.var=denom.var, time.res=time.res)
      names(q)<-geos
      return(q)
    }
    ) 
    names(ww)<- ages
    return(ww)
  }
  )
  names(all.glm.res)<-use.syndromes
  return(all.glm.res)
}

dashboardPlot<-function(all.glm.res){
  ds<-all.glm.res
  counties.to.test<- names(ds[[1]][[1]])
  ages.to.test<-names(ds[[1]])
  age.labels<-ages.to.test
  syndromes<-names(ds)
  dates<-as.Date(ds[[1]][[1]][[1]][['date']])
  n.times<-length(dates)
  last.date.format<-max(dates)
  last.date.format<-format(last.date.format,
                           "%b %d, %Y")
  

  server<-function(input, output){
    output$countyPlot = renderPlot({
      
    ili2.resid<- sapply(ds[[input$set.syndrome]], function(x) sapply(x,'[[','resid1'), simplify='array')
    dimnames(ili2.resid)[[2]]<-counties.to.test
    ili2.pred<- sapply(ds[[input$set.syndrome]], function(x) sapply(x,'[[','pred'), simplify='array')
    dimnames(ili2.pred)[[2]]<-counties.to.test
    ili2.pred.lcl<- sapply(ds[[input$set.syndrome]], function(x) sapply(x,'[[','lpi'), simplify='array')
    dimnames(ili2.pred.lcl)[[2]]<-counties.to.test
    ili2.pred.ucl<- sapply(ds[[input$set.syndrome]], function(x) sapply(x,'[[','upi'), simplify='array')
    dimnames(ili2.pred.ucl)[[2]]<-counties.to.test
    obs.ili<- sapply(ds[[input$set.syndrome]], function(x) sapply(x,'[[','y'), simplify='array')
    dimnames(obs.ili)[[2]]<-counties.to.test
    denom<- sapply(ds[[input$set.syndrome]], function(x) sapply(x,'[[','denom'), simplify='array')
    dimnames(denom)[[2]]<-counties.to.test
      plot.min<-which(input$display.dates==dates)
      dates.select<-plot.min:n.times
      par(mfrow=c(2,3), mar=c(3,2,1,1))
      for(i in ages.to.test){
        for( j in input$set.borough){
          if(input$set.prop=='Counts'){
            y=obs.ili[dates.select,j,i]
            pred<-ili2.pred[dates.select,j,i]
            pred.lcl<-ili2.pred.lcl[dates.select,j,i]
            pred.ucl<-ili2.pred.ucl[dates.select,j,i]
            if(input$set.axis==F){
              y.range<-c(0,max(c(ili2.pred.lcl[dates.select,j,i],ili2.pred.ucl[dates.select,j,i],ili2.pred[dates.select,j,i],obs.ili[dates.select,j,i]), na.rm=T))
            }else{
              y.range<-c(0,max(c(ili2.pred.lcl[dates.select,j,],ili2.pred.ucl[dates.select,j,],ili2.pred[dates.select,j,],obs.ili[dates.select,j,]), na.rm=T))
            }
          }else if (input$set.prop=='Proportion'){
            y=obs.ili[dates.select,j,i]/denom[dates.select,j,i]
            pred<-ili2.pred[dates.select,j,i]/denom[dates.select,j,i]
            pred.lcl<-ili2.pred.lcl[dates.select,j,i]/denom[dates.select,j,i]
            pred.ucl<-ili2.pred.ucl[dates.select,j,i]/denom[dates.select,j,i]
            
            if(input$set.axis==F){
              y.range<-c(0,max(y,na.rm=T))
            }else{
              y.range<-c(0, max(plot.prop[dates.select,j,], na.rm=T))
            }
          }else{
            y=obs.ili[dates.select,j,i]/ili2.pred[dates.select,j,i]
            pred<-obs.ili[dates.select,j,i]/ili2.pred[dates.select,j,i]
            pred.lcl<-obs.ili[dates.select,j,i]/ili2.pred.lcl[dates.select,j,i]
            pred.ucl<-obs.ili[dates.select,j,i]/ili2.pred.ucl[dates.select,j,i]
            if(input$set.axis==F){
              y.range<-range(y,na.rm=T)
              y.range[is.infinite(y.range)]<-10
            }else{
              y.range<-c(0.2, 4)
            }  
          }
          plot(dates[dates.select],y, type='n', bty='l', ylab='Fitted', main=paste(j, i), ylim=y.range)
          polygon(c(dates[dates.select],rev(dates[dates.select])), 
                  c(pred.lcl, rev(pred.ucl)), col=rgb(1,0,0,alpha=0.1), border=NA)
          lines(dates[dates.select],pred, type='l', col='red', lty=1, lwd=1.5 )
          lines(dates[dates.select],y, lwd=1.5)
          if(input$set.prop=='Observed/Expected'){
            abline(h=1, col='gray', lty=2)
          }
        }
      }
    }
    ,
    width = "auto", height = "auto")
  }
  
  ui<-
    fluidPage(
      titlePanel(paste0('Data through ', last.date.format)),
      sidebarLayout(
        sidebarPanel(
          selectInput("set.prop", "Proportion of ED visits or count:",
                      choice=c('Proportion','Counts','Observed/Expected'), selected ="Proportion" ),
          selectInput("set.borough", "Geographic unit:",
                      choice=counties.to.test, selected ="Citywide" ),
          selectInput("set.syndrome", "Syndrome:",
                      choice=syndromes, selected ="ili" ),
          checkboxInput("set.axis", "Uniform axis for all plots?:",
                        value =F ),
          sliderInput('display.dates', 'Earliest date to display', min=min(dates), step=7,max=dates[length(dates)-2], value=dates[length(dates)-round(length(dates)/5)]),
        ),
        mainPanel(
          plotOutput("countyPlot"),
          column(8, align = 'justify',
                 hr(),
                 span("The black line shows the observed number of ED visits per day in the indicated stratum, and the red lines denote the mean and 95% prediction intervals for a model adjusting for seasonality, influenza activity, and RSV activity"),
                 hr(),
                 span("This app and package were developed by The Public Health Modeling Unit and The Weinberger Lab at Yale School of Public Health. Contributors include Dan Weinberger, Alyssa Amick, Forrest Crawford, Kelsie Cassell, Marcus Rossi, Ernest Asare, Yu-Han Kao. Underlying analysis code can be found at https://github.com/weinbergerlab/ExcessILI"),
                 
          )
        )
      )
    )
  shinyApp(ui, server)
}

excessExtract<-function(ds, syndrome, extract.quantity){
  out.ds<- sapply(ds[[syndrome]], function(x) sapply(x,'[[',extract.quantity), simplify='array')
  return(out.ds)
}
