library(dplyr)
hts <- function(basis, sMat) {
  t(sMat %*% t(basis))
}
rmse <- function(a, b){
  rmses <- sqrt(colMeans((a-b)^2))
}

sMat = rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
             diag(rep(1, 4)))


for (file in c('ets.csv', 'ea.csv')){
  reconf = read.csv(paste0('scenario1_', file), row.names = 1)
  basef1 = read.csv('basef1.csv', row.names = 1)
  data1 = read.csv('data_scenario1.csv', row.names = 1)
  accs = data.frame()
  for (index in 1:100){
    acc = NULL
    observ = data1[data1$index==index & data1$t > 300, 1:4] %>% as.matrix() %>% unname()
    observ = hts(observ, sMat)
    # base forecast
    etsf = basef1[basef1$index==index & basef1$t > 300 & basef1$method=='ets', 1:7] %>% as.matrix() %>% unname()
    arimaf = basef1[basef1$index==index & basef1$t > 300 & basef1$method=='arima', 1:7] %>% as.matrix() %>% unname()
    acc = rbind(acc, rmse(etsf, observ), rmse(arimaf, observ))
    row_name = c('ets', 'arima')
    for (immu in 0:1){
      for (method in c('ols', 'wlss', 'wlsv', 'shrinkage')){
        fof = reconf[reconf$method==method & reconf$immu==immu & reconf$index==index, 1:7]
        acc = rbind(acc, rmse(fof, observ))
        row_name <- c(row_name, paste0(method, immu))
      }
    }
    acc = data.frame(acc)
    acc$method = row_name
    acc$index = index
    accs <- rbind(accs, acc)
  }
  write.csv(accs, paste0('acc_1_', file))
}


for (file in c('ets.csv', 'ea.csv')){
  reconf = read.csv(paste0('scenario2_', file), row.names = 1)
  basef1 = read.csv('basef2.csv', row.names = 1)
  data1 = read.csv('data_scenario2.csv', row.names = 1)
  accs = data.frame()
  for (index in 1:100){
    acc = NULL
    observ = data1[data1$index==index & data1$t > 300, 1:4] %>% as.matrix() %>% unname()
    observ = hts(observ, sMat)
    # base forecast
    etsf = basef1[basef1$index==index & basef1$t > 300 & basef1$method=='ets', 1:7] %>% as.matrix() %>% unname()
    arimaf = basef1[basef1$index==index & basef1$t > 300 & basef1$method=='arima', 1:7] %>% as.matrix() %>% unname()
    acc = rbind(acc, rmse(etsf, observ), rmse(arimaf, observ))
    row_name = c('ets', 'arima')
    for (immu in 0:1){
      for (method in c('ols', 'wlss', 'wlsv', 'shrinkage')){
        fof = reconf[reconf$index==index & reconf$method==method & reconf$immu==immu, 1:7]
        acc = rbind(acc, rmse(fof, observ))
        row_name <- c(row_name, paste0(method, immu))
      }
    }
    acc = data.frame(acc)
    acc$method = row_name
    acc$index = index
    accs <- rbind(accs, acc)
  }
  write.csv(accs, paste0('acc_2_', file))
}


panels = c('acc_1_ets.csv', 'acc_1_ea.csv', 'acc_2_ets.csv', 'acc_2_ea.csv')
for (i in 1:length(panels)){
  data = read.csv(panels[i], row.names = 1)
  data %>% group_by(method) %>% summarise(level1=mean(X1), 
                                          level2=mean(X2), 
                                          level3=mean(X3)) %>%
    t() %>% write.csv(paste0('panel', i, '.csv'))
}


