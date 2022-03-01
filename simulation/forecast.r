# utility functions
library(forecast)
source('../chf.r')

hts <- function(basis, sMat) {
  t(sMat %*% t(basis))
}

cal_basef <- function(series, method = 'arima'){
  series = hts(series, sMat)
  apply(series, 2, function(series){
    series = ts(series, frequency = 12)
    if (method == 'arima'){
      model = auto.arima(series)
    }
    if (method == 'ets'){
      model = ets(series)
    }
    c(fitted(model), forecast(model, h=24)$mean)
  })
}
reconcile <- function(series, basef, immu_set){
  series = hts(series, sMat)
  insample_f = basef[1:300,]
  insample_error = insample_f - series
  cov_mat = cov(insample_error)
  outsample_f = basef[301:324,]
  
  # ols 
  ols = forecast.reconcile(outsample_f, 
                     sMat, 
                     diag(rep(1, 7)),
                     immu_set = immu_set)
  
  # wlss
  weight_mat = diag(rowSums(sMat))
  wlss = forecast.reconcile(outsample_f, 
                            sMat, 
                            weight_mat,
                            immu_set = immu_set)
  # wlsv
  weight_mat = diag(diag(cov_mat))
  wlsv = forecast.reconcile(outsample_f, 
                           sMat, 
                           weight_mat,
                           immu_set = immu_set)
  # shrinkage
  rbind(data.frame(ols, method='ols', t=301:324),
        data.frame(wlss, method='wlss', t=301:324),
        data.frame(wlsv, method='wlsv', t=301:324))
}


# settings
sMat = rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
             diag(rep(1, 4)))


# scenario
data1 = read.csv('data_scenario1.csv', row.names = 1)
df = data.frame()
for (index in 1:100){
  print(index)
  series <- data1[data1$index==index,][1:300, 1:4]
  basef1 <- data.frame(cal_basef(series, method='arima'),
                      method='arima',
                      index=index,
                      t=1:324)
  basef2 <- data.frame(cal_basef(series, method='ets'),
                      method='ets',
                      index=index,
                      t=1:324)
  df <- rbind(df, basef1, basef2)
}
write.csv(df, 'basef1.csv')

# scenario 2
data2 = read.csv('data_scenario2.csv', row.names = 1)
df = data.frame()
for (index in 1:100){
  print(index)
  series <- data2[data2$index==index,][1:300, 1:4]
  basef1 <- data.frame(cal_basef(series, method='arima'),
                       method='arima',
                       index=index,
                       t=1:324)
  basef2 <- data.frame(cal_basef(series, method='ets'),
                       method='ets',
                       index=index,
                       t=1:324)
  df <- rbind(df, basef1, basef2)
}
write.csv(df, 'basef2.csv')

