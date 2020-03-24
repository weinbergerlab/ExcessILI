##Import auxillary data
#Pulls in google searches for RSV in CT can get 5 years of historical data at weekly reslution (switches to monthly for >5 years)
rsv.google.import<-function(state){
  rsv.down<-gtrends(keyword = 'RSV', geo=paste0("US-", state), time = "today+5-y",
                    gprop = c("web") , category = 0,  low_search_volume = FALSE)
  rsv<-rsv.down$interest_over_time[,c('date','hits')]
  rsv$date<-as.Date(rsv$date)
  mmwr.week.rsv<-MMWRweek(rsv$date)[,c('MMWRyear','MMWRweek')]
  rsv<-cbind.data.frame(rsv,mmwr.week.rsv)
  names(rsv)<-c('date','rsv.searches','MMWRyear','MMWRweek')
  return(rsv)
}

#Pull NREVSS testing data (% positive)
nrevss_flu_import<-function(){
  nrevvs.state<-who_nrevss(region = c("state"))
  clin<-nrevvs.state[['clinical_labs']]
  data(hhs_regions)
  cw.file<-hhs_regions
  clin2<-merge(clin, cw.file, by.x='region', by.y="state_or_territory")
  clin2<-clin2[,c('region','region_number','year','week','wk_date',"percent_positive")]
  names(clin2)[1:2]<-c('state','hhs_region')
  nrevvs_hhs<-who_nrevss(region = c("hhs"))
  clin.hhs<-nrevvs_hhs[['clinical_labs']]
  clin.hhs<-clin.hhs[,c('region','wk_date',"percent_positive")]
  clin.hhs$region<-as.numeric(gsub('Region ', '', clin.hhs$region))
  names(clin.hhs)<-c('hhs_region','wk_date',"hhs_percent_positive" )
  clin3<-merge(clin2, clin.hhs, by=c('hhs_region', 'wk_date'))
  clin3$percent_positive[is.na(clin3$percent_positive)]<-clin3$hhs_percent_positive[is.na(clin3$percent_positive)]
  clin3$state.abb<- state.abb[match(clin3$state,state.name)]
  names(clin3)<- c("hhs_region","wk_date", "state_name", "MMWRyear","MMWRweek" ,"flu_pct_pos","hhs_percent_positive", "state")
  clin3$flu_pct_pos<-as.numeric(as.character(clin3$flu_pct_pos))
  return(clin3)
}

combine_ds<-function(ds=combo2.sub, geo=geo, agevar=agevar, datevar=datevar){
ili.m<-melt(ds, id.vars=c(geo,agevar,datevar,'MMWRyear','MMWRweek'))
  form1<-as.formula(paste0())
  ili.a<-acast(ili.m, ddate+MMWRyear+MMWRweek ~geo~age~variable , fun.aggregate = sum )
  dimnames(ili.a)[[1]]<-substr(dimnames(ili.a)[[1]],1,10)
return(ili.a)
}


## Evaluate results after controlling for flu and RSV
glm.func<-function(ds, x.test, age.test, syndrome){
  date.string<-as.Date(dimnames(ds)[[1]])
  month<-month(date.string)
  epiyr<-year(date.string)
  epiyr[month<=6] <- epiyr[month<=6]-1
  epiyr.index<-epiyr-min(epiyr)+1
  weekN<-MMWRweek(date.string)[,'MMWRweek']
  day.of.year<-yday(date.string)
  day.of.week<-as.factor(weekdays(date.string))
  
  clean.array.citywide<-ds[,x.test,,]
  epiyr.index.f<-as.factor(epiyr.index)
  epiyr.index.f2<-as.factor(epiyr.index)
  
  y.age = t(clean.array.citywide[,age.test,syndrome])
  n.dates<-length(y.age)
  
  if(syndrome=='ili'){
    denom<-y.age[1,]/t(clean.array.citywide[,age.test,'ili.prop'])[1,]}else{
      denom<-y.age[1,]/t(clean.array.citywide[,age.test,'resp.prop']+0.01)[1,]
    }
  
  y.age.fit<-y.age[1,]
  y.age.fit[date.string>=as.Date("2020-03-01")] <- NA #extrapolate last 1 months
  
  sqrt.rsv =sqrt(clean.array.citywide[,age.test,'rsv.searches']) #same for all ages and boroughs
  sqrt.flu =sqrt(clean.array.citywide[,age.test,'flu_pct_pos']) #same for all ages and boroughs
  
  
  sqrt.flu<- na.locf(sqrt.flu)  #fill in missing observations for flu at end of TS with most recent observed values
  
  sqrt.rsv<- na.locf(sqrt.rsv)  #fill in missing observations for RSV at end of TS with most recent observed values
  
  t2<-1:length(y.age)
  # 
  sin1<-sin(2*pi*t2/365.25)
  cos1<-cos(2*pi*t2/365.25)
  sin2<-sin(2*pi*t2*2/365.25)
  cos2<-cos(2*pi*t2*2/365.25)
  sin3<-sin(2*pi*t2*3/365.25)
  cos3<-cos(2*pi*t2*3/365.25)
  
  
  ds.glm<-cbind.data.frame('day.of.year'=day.of.year,'y.age'=y.age[1,],'y.age.fit'=y.age.fit,sqrt.rsv, sqrt.flu, day.of.week, t2,epiyr.index.f, sin1, sin2, cos1, cos2, sin3, cos3)
  #ds.glm<-ds.glm[complete.cases(ds.glm),]
  form1<-as.formula(paste0('y.age.fit',"~ epiyr.index.f*sqrt.rsv +   #rsv effect varies by epiyr
                   epiyr.index.f*sqrt.flu + #flu effect, varies by epiyear
                   day.of.week+
                   sin1+cos1 +sin2+cos2+ sin3+cos3 "
  ))
  form2<-as.formula(paste0('y.age',"~ epiyr.index.f*sqrt.rsv +   #rsv effect varies by epiyr
                   epiyr.index.f*sqrt.flu + #flu effect, varies by epiyear
                   day.of.week+
                   sin1+cos1 +sin2+cos2+ sin3+cos3  "
  ))
  mod1<-glm(form1, data=ds.glm, family=poisson(link='log'))
  #500 samples total
  pred.coefs.reg.mean<- mvrnorm(n = 100, mu=coef(mod1), Sigma=vcov( mod1))
  mod.mat.pred<-model.matrix(form2, data=ds.glm, family='poisson')
  preds.stage1.regmean<- mod.mat.pred %*% t(pred.coefs.reg.mean) 
  preds.stage2<-rpois(n=length(preds.stage1.regmean)*5, exp(preds.stage1.regmean))
  preds.stage2<-matrix(preds.stage2, nrow=nrow(preds.stage1.regmean), ncol=ncol(preds.stage1.regmean)*5)
  preds.stage2.q<-t(apply(preds.stage2,1,quantile, probs=c(0.025,0.5, 0.975)))                   
  
  resid1<- log((ds.glm$y.age+0.5) /(preds.stage2.q[,'50%']+0.5))
  
  out.list<-list(y=y.age[1,], pred=preds.stage2.q[,'50%'], resid1=resid1, upi=preds.stage2.q[,'97.5%'], lpi=preds.stage2.q[,'2.5%'],'ili.prop'= t(clean.array.citywide[,age.test,'ili.prop']), 'resp.prop'= t(clean.array.citywide[,age.test,'resp.prop'])
  )
  return(out.list)
}


#SHINY dashboard functions
shiny.dashboard.func<-function(){
server<-function(input, output){
  output$countyPlot = renderPlot({
    ili.prop<-sapply(ds[[input$set.syndrome]], function(x) sapply(x,'[[','ili.prop'), simplify='array')
    resp.prop<-sapply(ds[[input$set.syndrome]], function(x) sapply(x,'[[','resp.prop'), simplify='array')
    dimnames(ili.prop)[[2]]<-counties.to.test
    dimnames(resp.prop)[[2]]<-counties.to.test
    if(input$set.syndrome=='ili'){
      plot.prop<-ili.prop
    }else{
      plot.prop<-resp.prop
    }
    
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
    #dates<-as.Date(dimnames(ili.a)[[1]])
    age.labels = c("Ages 0-4 years", "Ages 5-17 years", "Ages 18-64 years", "Ages 65+ years", "All age groups")
    
    plot.min<-which(input$display.dates==dates)
    dates.select<-plot.min:n.times
    par(mfrow=c(2,3), mar=c(3,2,1,1))
    for(i in c('1','2','3','4','5')){
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
        }else if (input$set.prop=='Counts/100,000 people'){
          y=obs.ili[dates.select,j,i]/pop3[dates.select,j,i]*100000
          pred<-ili2.pred[dates.select,j,i]/pop3[dates.select,j,i]*100000
          pred.lcl<-ili2.pred.lcl[dates.select,j,i]/pop3[dates.select,j,i]*100000
          pred.ucl<-ili2.pred.ucl[dates.select,j,i]/pop3[dates.select,j,i]*100000
          if(input$set.axis==F){
            y.range<-c(0,max(c(pred.lcl,pred.ucl,pred,y), na.rm=T))
          }else{
            y.range<-c(0,max(c(ili2.pred.lcl[dates.select,j,]/pop3[dates.select,j,]*100000,
                               ili2.pred.ucl[dates.select,j,]/pop3[dates.select,j,]*100000,
                               ili2.pred[dates.select,j,]/pop3[dates.select,j,]*100000,
                               obs.ili[dates.select,j,]/pop3[dates.select,j,]*100000
            ), na.rm=T))
          }
          
        }else if (input$set.prop=='Proportion'){
          y=plot.prop[dates.select,j,i]
          denom<-obs.ili[dates.select,j,i]/y
          pred<-ili2.pred[dates.select,j,i]/denom
          pred.lcl<-ili2.pred.lcl[dates.select,j,i]/denom
          pred.ucl<-ili2.pred.ucl[dates.select,j,i]/denom
          
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
        plot(dates[dates.select],y, type='n', bty='l', ylab='Fitted', main=paste(j, age.labels[as.numeric(i)]), ylim=y.range)
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
    titlePanel(paste0('NYC ED syndromic surveillance through ', last.date.format)),
    span("CAUTION: Syndromic surveillance data can be hard to interpret. Any increases above expected could be due to changes in healthcare seeking behavior (people might be more likely to go to the ED now with less severe symptoms because they are aware of the COVID-19 epidemic), or it could be due to actual viral illness, or a combination. For a deep dive of the data produced by NYC Department of Health and Mental Hygiene see https://www1.nyc.gov/assets/doh/downloads/pdf/hcp/weekly-surveillance03072020.pdf ."),
    sidebarLayout(
      sidebarPanel(
        selectInput("set.prop", "Proportion of ED visits or count:",
                    choice=c('Proportion','Counts','Counts/100,000 people','Observed/Expected'), selected ="Proportion" ),
        selectInput("set.borough", "Borough:",
                    choice=counties.to.test, selected ="Citywide" ),
        selectInput("set.syndrome", "Syndrome:",
                    choice=syndromes, selected ="ili" ),
        checkboxInput("set.axis", "Uniform axis for all plots?:",
                      value =F ),
        sliderInput('display.dates', 'Earliest date to display', min=min(dates), max=max(dates)-30, value=max(dates)-180),
      ),
      mainPanel(
        plotOutput("countyPlot"),
        column(8, align = 'justify',
               hr(),
               span("The black line shows the observed number of ED visits per day in the indicated stratum, and the red lines denote the mean and 95% prediction intervals for a model adjusting for seasonality, influenza activity, and RSV activity"),
               hr(),
               span("These plots summarize the NYC syndromic surveillance data, which were downloaded from the Epiquery website of the NYC Department of Health and Mental Hygiene. The models and plots were done by Dr. Dan Weinberger from the Public Health Modeling Unit and Department of Epidemiology of Microbial Diseases at Yale School of Public Health, with assistance from Alyssa Amick, Forrest Crawford, Kelsie Cassell, Marcus Rossi, Ernest Asare, Yu-Han Kao. Underlying analysis code can be found at https://github.com/weinbergerlab/covid19_syndromic"),
               
        )
      )
    )
  )
}

