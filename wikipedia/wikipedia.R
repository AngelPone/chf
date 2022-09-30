library(forecast)
library(foreach)
library(dplyr)
cluster_cores = 7
cl <- parallel::makeCluster(cluster_cores)
doParallel::registerDoParallel(cl)

get_forecasts <- function(x, index){
  library(forecast)
  x <- ts(x, frequency = 7)
  if (sum(x[307:366] == 0) > 36){
    m <- ses(x, h = 28)
    method <- "ses"
  } else {
    m <- forecast(ets(x), h=28)
    method <- "ets"
  }
  fcasts <- data.frame(matrix(c(m$fitted, m$mean), nrow = 1))
  colnames(fcasts) <- paste0('d_', 1:dim(fcasts)[2])
  fcasts$index <- index
  fcasts$method <- method
  fcasts
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

transform.sMat <- function(sMat, basis_set){
  m <- dim(sMat)[2]
  if (length(basis_set) != m){
    stop(simpleError(sprintf('length of basis set should be %d', m)))
  }
  S1 <- sMat[basis_set,]
  S2 <- sMat[-basis_set,]
  transitionMat <- solve(S1, diag(rep(1, m)))
  rbind(S2 %*% transitionMat, diag(rep(1, m)))
}

forecast.reconcile <- function(base_forecasts, 
                               sMat,
                               weighting_matrix,
                               immu_set=NULL,
                               nonnegative=FALSE){
  m <- dim(sMat)[2]
  n <- dim(sMat)[1]
  k <- length(immu_set)
  if (length(immu_set) == 0){
    weighting_matrix = solve(weighting_matrix)
    solution <- matrix(0, dim(base_forecasts)[1], m)
    if (nonnegative){
      for (i in 1:dim(base_forecasts)[1]){
        Dmat <- t(sMat) %*% weighting_matrix %*% sMat
        dvec <- t(sMat) %*% weighting_matrix %*% base_forecasts[i,]
        Amat <- diag(rep(1, m))
        solution[i,] <- quadprog::solve.QP(Dmat, dvec, Amat)$solution
      }
      reconciled_y <- sMat %*% t(solution)
    }else{
      reconciled_y <- sMat %*% solve(t(sMat) %*% weighting_matrix %*% sMat) %*% t(sMat) %*% weighting_matrix %*% t(base_forecasts)
    }
    return(t(reconciled_y))
  }
  
  # construct new basis time series
  if (k > m) {
    stop(simpleError(sprintf('length of basis set can not be bigger than %d', m)))
  }
  ## select mutable series
  immutable_basis <- sort(immu_set)
  candidate_basis <- setdiff((n-m+1):n, immu_set)
  determined <- 1:(n-m)
  mutable_basis <- candidate_basis
  if (any(immutable_basis < n-m+1)){
    i <- max(which(immutable_basis < n-m+1))
    while (i > 0) {
      corresponding_leaves <- which(sMat[immutable_basis[i], ] != 0) + n - m
      free_leaves <- setdiff(corresponding_leaves, c(immutable_basis, determined))
      if (length(free_leaves) == 0){
        if (all(corresponding_leaves %in% immutable_basis)){
          warning(paste0("all children of ", immutable_basis[i], "th series are immutable, it is removed from the condition."))
          k <- k - 1
          immutable_basis <- immutable_basis[immutable_basis != immutable_basis[i]]
          i <- i - 1
          next
        } else{
          stop(simpleError("can not describe the hierarchy."))
        }
      }
      determined <- determined[determined != immutable_basis[i]]
      determined <- c(determined, free_leaves[1])
      mutable_basis <- mutable_basis[mutable_basis != free_leaves[1]]
      i <- i - 1
    }
  }
  new_basis <- c(sort(mutable_basis), immutable_basis)
  sMat <- transform.sMat(sMat, new_basis)
  S1 <- sMat[1:(n-k),,drop=FALSE][,1:(m-k),drop=FALSE]
  S2 <- sMat[1:(n-m),,drop=FALSE][,(m-k+1):m,drop=FALSE]
  determined <- setdiff(1:n, new_basis)
  mutable_series <- c(determined, mutable_basis)
  mutable_weight <- solve(weighting_matrix)[mutable_series,,drop=FALSE][,mutable_series,drop=FALSE]
  mutable_base <- cbind(base_forecasts[,determined,drop=FALSE] - t(S2 %*% t(base_forecasts[,immutable_basis,drop=FALSE])),
                        base_forecasts[,sort(mutable_basis),drop=FALSE])
  reconciled_mutable <- solve(t(S1) %*% mutable_weight %*% S1) %*% t(S1) %*% mutable_weight %*% t(mutable_base)
  reconciled_y <- t(sMat %*% rbind(reconciled_mutable, t(base_forecasts[,immutable_basis,drop=FALSE])))
  mutable_weight <- mutable_weight / max(diag(mutable_weight))
  if (nonnegative){
    for (i in 1:dim(mutable_base)[1]){
      Dmat <- t(S1) %*% mutable_weight %*% S1
      dvec <- as.vector(t(mutable_base[i,]) %*% mutable_weight %*% S1)
      Amat <- diag(rep(1, dim(S1)[2]))
      bvec <- rep(0, dim(S1)[2])
      sol <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec)$solution)
      if (is(sol, "try-error")){
        warning(paste0("unsolvable at row ", rownames(basef)[i], " use unconstrained solution!"))
      }else{
        reconciled_y[i,] <- as.vector(sMat %*% c(sol, base_forecasts[i,immutable_basis,drop=FALSE]))
      }
    }
  }
  new_index <- c(determined, new_basis)
  reconciled_y[,order(new_index)]
}

reconcile_f <- function(basef, immu_set, nn=FALSE){
  
  stopifnot(all(dim(basef) == c(394, n)))
  basef_insample <- basef[1:366,]
  basef_outsample <- basef[367:394,]
  
  ols_weight <- diag(n)
  ols <- forecast.reconcile(basef_outsample, smat, ols_weight, immu_set = immu_set, nonnegative = nn)
  wlss_weight <- diag(as.vector(smat %*% matrix(1, nrow = m)))
  wlss <- forecast.reconcile(basef_outsample, smat, wlss_weight, immu_set = immu_set, nonnegative = nn)
  
  error_cov <- var(basef_insample - train)
  wlsv_weight <- diag(diag(error_cov))
  wlsv <- forecast.reconcile(basef_outsample, smat, wlsv_weight, immu_set = immu_set, nonnegative = nn)
  
  lamb <- lambda_estimate(basef_insample - train)
  shrink_weight <- (1-lamb) * error_cov + lamb * wlsv_weight
  shrink <- forecast.reconcile(basef_outsample, smat, shrink_weight, immu_set = immu_set, nonnegative = nn)
  
  list(basef=basef_outsample, ols=ols, wlss=wlss, wlsv=wlsv, shrink=shrink)
}

evaluate_f <- function(fcasts){
  level <- strsplit(colnames(wikipedia), '_') %>% sapply(function(x){x[1]})
  
  res <- sapply(fcasts, function(x){
    sqrt(colMeans((x - test)^2)) %>%
      split(level) %>%
      sapply(mean)
  }) %>% 
    as.data.frame() %>%
    mutate(level = row.names(.))
  row.names(res) <- NULL
  res
}





# 
wikipedia <- read.csv("wikipedia/wikimedia2016_hierarchy.csv")
train <- wikipedia[1:366,] %>% as.matrix() %>%unname()
test <- wikipedia[367:394,] %>% as.matrix()  %>% unname()
smat <- unname(as.matrix(read.csv('wikipedia/smat.csv')))
n <- dim(smat)[1]
m <- dim(smat)[2]

basef <- foreach(i=1:dim(wikipedia)[2], .combine = rbind) %dopar%{
  get_forecasts(as.numeric(wikipedia[1:(394-28), i]), index=i-1)
}


immutable_set <- basef %>% filter(method == "ses") %>%
  pull(index)
immutable_set <- immutable_set + 1
rec_f <- basef %>%
  arrange(index) %>%
  select(starts_with('d_')) %>%
  as.matrix() %>%
  unname() %>%
  t() %>%
  reconcile_f(c(1, immutable_set))

rec_f_U <- basef %>%
  arrange(index) %>%
  select(starts_with('d_')) %>%
  as.matrix() %>%
  unname() %>%
  t() %>%
  reconcile_f(immu_set = NULL)

rec_f_nn <- basef %>%
  arrange(index) %>%
  select(starts_with('d_')) %>%
  as.matrix() %>%
  unname() %>%
  t() %>%
  reconcile_f(c(1, immutable_set), nn=TRUE)

rec_f_U_nn <- basef %>%
  arrange(index) %>%
  select(starts_with('d_')) %>%
  as.matrix() %>%
  unname() %>%
  t() %>%
  reconcile_f(immu_set = NULL, nn=TRUE)

# accuracy

acc <- evaluate_f(rec_f) %>% as.data.frame() %>%
  mutate(constraints = "C")

acc <- evaluate_f(rec_f_U) %>% as.data.frame() %>%
  mutate(constraints = "U") %>%
  rbind(acc)

acc_nn <- evaluate_f(rec_f_nn) %>% as.data.frame() %>%
  mutate(constraints = "C")

acc_nn <- evaluate_f(rec_f_U_nn) %>% as.data.frame() %>%
  mutate(constraints = "U") %>%
  rbind(acc_nn)


acc[,1:5] <- round(acc[,1:5])
acc_nn[, 1:5] <- round(acc_nn[,1:5])

# write results to CSV

level_order = c("total", "group", "network", "language", 
                "language.group", "language.network",
                "access", "access.group", "access.network",
                "access.language", "access.language.group", 
                "access.language.network")

acc %>% 
  tidyr::pivot_wider(
  id_cols = "level",
  names_from = c("constraints"), 
  values_from = c("basef", "ols", "wlss", "wlsv", "shrink")) %>%
  arrange(factor(level, levels = level_order)) %>%
  write.csv('wikipedia/accuracy.csv')

acc_nn %>%
  tidyr::pivot_wider(
  id_cols = "level",
  names_from = c("constraints"), 
  values_from = c("basef", "ols", "wlss", "wlsv", "shrink")) %>%
  arrange(factor(level, levels = level_order)) %>%
  write.csv('wikipedia/accuracy_nn.csv')

