# utility functions
cumsum_matrix <- function(x){ apply(x, 2, cumsum) }
brodcast_sum <- function(x, y){
  t(apply(x, 1, function(row){row + y}))
}

# simulate time series given parameters in scenario1
simulate_ts <- function(n, sigma2_varrho, sigma2_varsigma, sigma2_omega, eta_mat, m=1){
  library(MASS)
  library(MTS)
  nseries <- dim(sigma2_omega)[1]
  
  # level
  initial_level <- rnorm(nseries)
  varrho <- mvrnorm(n-1, rep(0, nseries), sigma2_varrho)
  level <- brodcast_sum(cumsum_matrix(varrho), initial_level)
  level <- rbind(initial_level, level)
  
  # trend
  initial_trend <- rnorm(nseries)
  varsigma <- mvrnorm(n-1, rep(0, nseries), sigma2_varsigma)
  trend <- rbind(initial_trend, 
                 brodcast_sum(cumsum_matrix(varsigma), initial_trend))
  
  # seasonal
  initial_seasonal <- mvrnorm(m-1, rep(0, nseries), diag(rep(1, nseries)))
  omega <- mvrnorm(n-m+1, rep(0, nseries), sigma2_omega)
  seasons <- initial_seasonal
  for (i in m:n) {
    seasons <- rbind(seasons, -apply(seasons[(i-m+1):(i-1),], 2, sum)+omega[i-m+1,])
  }
  # eta
  eta <- VARMAsim(n, sample(c(0, 1), 1, replace = TRUE), sample(c(0, 1), 1, replace = TRUE),
                  phi = diag(runif(4, 0.5, 0.7)), theta = diag(runif(4, 0.5, 0.7)),
                  sigma = eta_mat)$series
  bts <- level + trend + seasons + eta
  bts
}

# simulate time series given parameters in scenario2
simulate_ts2 <- function(n, sigma2_varrho, sigma2_varsigma, sigma2_omega, eta_mat, error1, error2, m=1){
  library(MASS)
  library(MTS)
  nseries <- dim(sigma2_omega)[1]
  
  # level
  initial_level <- rnorm(nseries)
  varrho <- mvrnorm(n-1, rep(0, nseries), sigma2_varrho)
  level <- brodcast_sum(cumsum_matrix(varrho), initial_level)
  level <- rbind(initial_level, level)
  
  # trend
  initial_trend <- rnorm(nseries)
  varsigma <- mvrnorm(n-1, rep(0, nseries), sigma2_varsigma)
  trend <- rbind(initial_trend, 
                 brodcast_sum(cumsum_matrix(varsigma), initial_trend))
  
  # seasonal
  initial_seasonal <- mvrnorm(m-1, rep(0, nseries), diag(rep(1, nseries)))
  omega <- mvrnorm(n-m+1, rep(0, nseries), sigma2_omega)
  seasons <- initial_seasonal
  for (i in m:n) {
    seasons <- rbind(seasons, -apply(seasons[(i-m+1):(i-1),], 2, sum)+omega[i-m+1,])
  }
  # eta
  eta <- VARMAsim(n, sample(c(0, 1), 1, replace = TRUE), sample(c(0, 1), 1, replace = TRUE),
                  phi = diag(runif(4, 0.5, 0.7)), theta = diag(runif(4, 0.5, 0.7)),
                  sigma = eta_mat)$series
  bts <- level + trend + seasons + eta
  e1 <- rnorm(n, 0, sqrt(error1))
  e2 <- rnorm(n, 0, sqrt(error2))
  bts[, 1] = bts[, 1] - e1 - 0.5*e2
  bts[, 2] = bts[, 2] + e1 - 0.5*e2
  bts[, 3] = bts[, 3] - e1 + 0.5*e2
  bts[, 4] = bts[, 4] + e1 + 0.5*e2
  bts
}


# set parameters
set.seed(42)
n = 324
repeat_n = 100
freq = 12
sigma2_varrho = diag(rep(2, 4))
sigma2_varsigma = diag(rep(0.007, 4))
sigma2_varomega = diag(rep(7, 4))
eta_mat = diag(rep(3, 4))
eta_mat[1, 2] = eta_mat[2, 1] = -2
eta_mat[3, 4] = eta_mat[4, 3] = -1
error1 = 10
error2 = 9

# scenario1
simulated_set = data.frame()
for (i in 1:repeat_n) {
  ts = simulate_ts(n, 
                   sigma2_varrho, 
                   sigma2_varsigma, 
                   sigma2_varomega, 
                   eta_mat,
                   freq)
  ts = data.frame(ts, t=1:n, index=i, row.names = NULL)
  simulated_set = rbind(simulated_set, ts)
}

write.csv(simulated_set, file = 'data_scenario1.csv')

# scenario2
simulated_set = data.frame()
for (i in 1:repeat_n) {
  ts = simulate_ts2(n, 
                   sigma2_varrho, 
                   sigma2_varsigma, 
                   sigma2_varomega, 
                   eta_mat,
                   error1,
                   error2,
                   freq)
  ts = data.frame(ts, t=1:n, index=i, row.names = NULL)
  simulated_set = rbind(simulated_set, ts)
}

write.csv(simulated_set, file = 'data_scenario2.csv')
