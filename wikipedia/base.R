wikipedia <- read.csv('wikimedia2016_hierarchy.csv')

library(forecast)
library(foreach)
cluster_cores = 7
cl <- parallel::makeCluster(cluster_cores)
doParallel::registerDoParallel(cl)

get_forecasts <- function(x, index){
  library(forecast)
  x <- ts(x, frequency = 7)
  m_arima <- forecast(auto.arima(x), h=28)
  m_ets <- forecast(ets(x), h=28)
  fcasts <- data.frame(rbind(c(m_arima$fitted, m_arima$mean), c(m_ets$fitted, m_ets$mean)))
  colnames(fcasts) <- paste0('d_', 1:dim(fcasts)[2])
  fcasts$index <- index
  fcasts$method <- c('arima', 'ets')
  fcasts
}
res <- foreach(i=1:dim(wikipedia)[2], .combine = rbind) %dopar%{
  get_forecasts(as.numeric(wikipedia[1:(394-28), i]), index=i-1)
}

write.csv(res, 'base.csv', row.names=FALSE)
