# source("./R/functions.R")
# source("./R/aux_functions.R")

# Open Rproj, which sets Working directory automatically)
combo1 <- readRDS("./Data/NYC_resp_ili_combo.rds")

combo1$denom1 <- combo1$ili  / (combo1$ili.prop  + 0.01)
combo1$denom2 <- combo1$resp / (combo1$resp.prop + 0.01)

combo1$total.visits <- round((combo1$denom1 + combo1$denom2)/2) + 1

combo1$state <- "NY"
combo1 <- combo1[combo1$borough %in% c("Citywide", "Bronx", "Manhattan") &
                 combo1$agec %in% c("1", "2", "3", "4", "5"), ]

combo1$agec <- as.character(combo1$agec)
combo1$agec[combo1$agec == "1"] <- "i. u5y"
combo1$agec[combo1$agec == "2"] <- "ii. 5-17y"
combo1$agec[combo1$agec == "3"] <- "iii. 18-64y"
combo1$agec[combo1$agec == "4"] <- "iv. 65+y"
combo1$agec[combo1$agec == "5"] <- "v. All ages"

excess_cases1 <-
  excessCases(ds = combo1,
              sub.statevar  = "borough",
              datevar       = "ddate", 
              agevar        = "agec",
              statevar      = "state",
              denom.var     = "total.visits",
              use.syndromes = c("ili", "resp"),
              adj.flu       = F,
              flu.import    = F)

dashboardPlot(excess_cases1)

ili2.resid <- sapply(excess_cases1[['ili']],
                     function(x) sapply(x, "[[", "resid1"), 
                     simplify = "array")

unexplained.cases <-
  excessExtract(ds = excess_cases1,
                syndrome = "ili",
                extract.quantity = "unexplained.cases")

excess.rr <-
  excessExtract(ds = excess_cases1,
                syndrome = "ili",
                extract.quantity = "resid1")


par(mfrow=c(1,1))
matplot(exp(excess.rr[,,1]), type='l')
matplot(unexplained.cases[,-2,1], type='l')

#Try out ts format package
n.obs<-10000
set.seed(123)
sim1<-as.data.frame(matrix(NA,nrow=n.obs, ncol=5))
names(sim1)<-c('state','date','agegrp','ili','resp')
sim1$state<- c(rep('CT', times=n.obs*0.3), rep("NY", times=n.obs*0.7))
sim1$agegrp<-sample(1:5, n.obs, replace=T)
sim1$date<-sample(seq.Date(from=as.Date('2019-01-01'),by='day', length.out=500), 1000, replace=T)
sim1$ili<-rbinom(n=n.obs, size=1, prob=0.05)
sim1$resp<-rbinom(n=n.obs, size=1, prob=0.1)
ts1<-ts_format(line.list=sim1, datevar='date',agevar='agegrp',statevar='state', syndromes=c('ili','resp'))

