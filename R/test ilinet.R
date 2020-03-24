library(lubridate)
library(MMWRweek)
library(cdcfluview)
library(gtrendsR)
library(shiny)
library(pbapply)
library(zoo)
library(MASS)
library(tidyr)
library(reshape2)

ili.data<-ilinet(region = c("state"))
ili.data$state<- state.abb[match(ili.data$region,state.name)]
ili.data<-ili.data[,c('state','week_start', 'ilitotal','total_patients')]
ili.data<-ili.data[ili.data$state %in% c('CA','NY','WA','NJ','CT'),]

source('./R/functions.R')
source('./R/aux_functions.R')

excess_cases1<-excessCases(ds=ili.data, 
                           datevar='week_start',
                           agevar='agec',
                           statevar='state',
                           denom.var="total_patients",
                           use.syndromes=c("ilitotal"), 
                           adj.flu=F, 
                           flu.import=F,
                           time.res='week')
dashboardPlot(excess_cases1)
