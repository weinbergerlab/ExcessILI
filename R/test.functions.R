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

#Open Rproj, which sets Working directory automatically)
combo1<-readRDS( './Data/NYC_resp_ili_combo.rds')
combo1$denom1<-combo1$ili/(combo1$ili.prop+0.01)
combo1$denom2<-combo1$resp/(combo1$resp.prop+0.01)
combo1$total.visits<-round((combo1$denom1+combo1$denom2)/2)+1
combo1$state<-'NY'
combo1<-combo1[combo1$borough %in% c('Citywide','Bronx','Manhattan') & combo1$agec %in% c('1','2','3','4','5'),]

source('./R/functions.R')
source('./R/aux_functions.R')

excess_cases1<-excessCases(ds=combo1, 
                           geovar='borough',
                           datevar='ddate',
                           agevar='agec',
                           statevar='state',
                           denom.var='total.visits',
                           use.syndromes=c('ili','resp') )
