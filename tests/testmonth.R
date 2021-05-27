library(ExcessILI)

d1 <- read.csv('./tests/Underlying Cause of Death, 1999-2019.csv', header=TRUE)

 d1$date <- as.Date(paste(substr(d1$Month.Code,1,4),substr(d1$Month.Code,6,7),'01', sep='-' ))

 d1 <- d1[,c('date','State','Deaths')]
 d1 <- d1[order(d1$State, d1$date),]
 d1 <- d1[!is.na(d1$date),]

d1$one <- 1
d1 <- d1[d1$State %in% c('New York','California','Texas'),]
run1 <- excessCases(d1, statevar='State',datevar = 'date', time.res='month',
                    adj.month.days=T, 
                    extrapolation.date = as.Date("2019-04-01"),
                    sum.dates = as.Date("2019-04-01"),
                    use.syndromes='Deaths', denom.var = 'one',extend.epiyear = T)

run2 <- excessCases(d1, statevar='State',datevar = 'date', time.res='month',
                    adj.month.days=F, 
                    extrapolation.date = as.Date("2019-03-01"),
                    sum.dates = as.Date("2019-03-01"),
                    use.syndromes='Deaths', denom.var = 'one',extend.epiyear = T)


str(run1)

plot(run1$Deaths$`1`$California$date,run1$Deaths$`1`$California$y)
points(run1$Deaths$`1`$California$date,run1$Deaths$`1`$California$pred, type='l')


plot(run2$Deaths$`1`$California$date,run2$Deaths$`1`$California$y)
points(run2$Deaths$`1`$California$date,run2$Deaths$`1`$California$pred, type='l')
