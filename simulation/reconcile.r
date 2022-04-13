source('chf.r')
library(dplyr)
hts <- function(basis, sMat) {
  t(sMat %*% t(basis))
}

lambda_estimate <- function(error){
  timeT = dim(error)[1]
  covm = t(error) %*% error/timeT
  xs = apply(error, 1, function(x){x / sqrt(diag(covm))}) %>% t()
  corm = t(xs) %*% xs / timeT
  diag(corm) = 0
  d = sum(corm^2)
  xs2 = xs^2
  v = 1/(timeT*(timeT-1))*(t(xs2) %*% xs2 - 1/timeT*(t(xs) %*% xs)^2)
  diag(v) = 0
  max(min(c(sum(v)/d, 1)), 0)
}

reconcile <- function(series, basef, immu_set){
  series = hts(series, sMat)
  insample_f = basef[1:300,]
  insample_error = insample_f - series
  cov_mat = cov(insample_error)
  outsample_f = basef[301:324,]
  write.csv(insample_error, 'test.csv')
  # conventional
  ols_conv = forecast.reconcile(outsample_f, 
                            sMat, 
                            diag(rep(1, 7)),
                            immu_set = NULL)
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
  wlss_conv = forecast.reconcile(outsample_f, 
                            sMat, 
                            weight_mat,
                            immu_set = NULL)
  # wlsv
  weight_mat = diag(diag(cov_mat))
  wlsv = forecast.reconcile(outsample_f, 
                            sMat, 
                            weight_mat,
                            immu_set = immu_set)
  wlsv_conv = forecast.reconcile(outsample_f, 
                            sMat, 
                            weight_mat,
                            immu_set = NULL)
  # shrinkage
  lamb = lambda_estimate(insample_error)
  weight_mat = (1-lamb) * cov_mat + lamb * weight_mat
  shrinkage = forecast.reconcile(outsample_f, 
                                 sMat, 
                                 weight_mat,
                                 immu_set = immu_set)
  shrinkage_conv = forecast.reconcile(outsample_f, 
                            sMat, 
                            weight_mat,
                            immu_set = NULL)
  
  rbind(data.frame(ols, method='ols', t=301:324, immu=1),
        data.frame(wlss, method='wlss', t=301:324, immu=1),
        data.frame(wlsv, method='wlsv', t=301:324, immu=1),
        data.frame(shrinkage, method='shrinkage', t=301:324, immu=1),
        data.frame(ols_conv, method='ols', t=301:324, immu=0),
        data.frame(wlss_conv, method='wlss', t=301:324, immu=0),
        data.frame(wlsv_conv, method='wlsv', t=301:324, immu=0),
        data.frame(shrinkage_conv, method='shrinkage', t=301:324, immu=0))
}

# load data
data1 = read.csv('data_scenario1.csv', row.names = 1)
basef1 = read.csv('basef1.csv', row.names = 1)

recons1 = data.frame()
for (index in 1:100){
  foo1 = as.matrix(data1[data1$t<=300 & data1$index == index, 1:4]) %>%
    unname()
  foof1 = as.matrix(basef1[basef1$index==index & basef1$method=='ets', 1:7]) %>% 
    unname()
  recon1 = reconcile(foo1, foof1, c(1))
  recon1$index = index
  recons1 = rbind(recons1, recon1)
}
write.csv(recons1, 'scenario1_ets.csv')


recons1 = data.frame()
for (index in 1:100){
  foo1 = as.matrix(data1[data1$t<=300 & data1$index == index, 1:4]) %>%
    unname()
  foof11 = as.matrix(basef1[basef1$index==index & basef1$method=='ets', 1:1, drop=FALSE]) %>% 
    unname()
  foof12 = as.matrix(basef1[basef1$index==index & basef1$method=='arima', 2:7, drop=FALSE]) %>% 
    unname()
  foof1 = cbind(foof11, foof12)
  recon1 = reconcile(foo1, foof1, c(1))
  recon1$index = index
  recons1 = rbind(recons1, recon1)
}
write.csv(recons1, 'scenario1_ea.csv')



data1 = read.csv('data_scenario2.csv', row.names = 1)
basef1 = read.csv('basef2.csv', row.names = 1)

recons1 = data.frame()
for (index in 1:100){
  foo1 = as.matrix(data1[data1$t<=300 & data1$index == index, 1:4]) %>%
    unname()
  foof1 = as.matrix(basef1[basef1$index==index & basef1$method=='ets', 1:7]) %>% 
    unname()
  recon1 = reconcile(foo1, foof1, c(1))
  recon1$index = index
  recons1 = rbind(recons1, recon1)
}
write.csv(recons1, 'scenario2_ets.csv')



recons1 = data.frame()
for (index in 1:100){
  foo1 = as.matrix(data1[data1$t<=300 & data1$index == index, 1:4]) %>%
    unname()
  foof11 = as.matrix(basef1[basef1$index==index & basef1$method=='ets', 1:1, drop=FALSE]) %>% 
    unname()
  foof12 = as.matrix(basef1[basef1$index==index & basef1$method=='arima', 2:7, drop=FALSE]) %>% 
    unname()
  foof1 = cbind(foof11, foof12)
  recon1 = reconcile(foo1, foof1, c(1))
  recon1$index = index
  recons1 = rbind(recons1, recon1)
}
write.csv(recons1, 'scenario2_ea.csv')






