
# produce summing matrix of new basis time series
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


# 
forecast.reconcile <- function(base_forecasts, 
                               sMat,
                               weighting_matrix,
                               immu_set=NULL){
  m <- dim(sMat)[2]
  n <- dim(sMat)[1]
  k <- length(immu_set)
  if (length(immu_set) == 0){
    weighting_matrix = solve(weighting_matrix)
    reconciled_y = sMat %*% solve(t(sMat) %*% weighting_matrix %*% sMat) %*% t(sMat) %*% weighting_matrix %*% t(base_forecasts)
    return(t(reconciled_y))
  }

  # construct new basis time series
  if (k > m) {
    stop(simpleError(sprintf('length of basis set can not be bigger than %d', m)))
  }
  ## select mutable series
  candidate_basis <- setdiff((n-m+1):n, immu_set)
  mutable_basis <- c()
  immutable_basis <- sort(immu_set)
  determined <- c()
  i <- max(which(immutable_basis < n-m+1))
  if (i > 0){
    while (length(mutable_basis) != m-k) {
      corresponding_leaves <- which(sMat[immutable_basis[i], ] == 1) + n - m
      free_leaves <- setdiff(corresponding_leaves, c(immutable_basis, mutable_basis, determined))
      if (length(free_leaves) == 0) stop(simpleError('the immu_set can not be used to describe the hierarchy'))
      if (length(free_leaves) == 1) {
        candidate_basis <- candidate_basis[candidate_basis != free_leaves[1]]
      } else{
        determined <- c(determined, free_leaves[1])
        mutable_basis <- c(mutable_basis, free_leaves[2:length(free_leaves)])
        candidate_basis <- candidate_basis[!(candidate_basis %in% free_leaves)]
      }
      i <- i - 1
    }
  } else {
    mutable_basis <- setdiff((n-m+1):n, immu_set)
  }
  
  new_basis <- c(sort(mutable_basis), immutable_basis)
  sMat <- transform.sMat(sMat, new_basis)
  S1 <- sMat[1:(n-k),,drop=FALSE][,1:(m-k),drop=FALSE]
  S2 <- sMat[1:(n-m),,drop=FALSE][,(m-k+1):m,drop=FALSE]
  determined <- setdiff(1:n, new_basis)
  mutable_series <- c(determined, mutable_basis)
  mutable_weight <- solve(weighting_matrix[mutable_series,,drop=FALSE][,mutable_series,drop=FALSE])
  mutable_base <- cbind(base_forecasts[,determined,drop=FALSE] - t(S2 %*% t(base_forecasts[,immutable_basis,drop=FALSE])),
                        base_forecasts[,sort(mutable_basis),drop=FALSE])
  reconciled_mutable <- solve(t(S1) %*% mutable_weight %*% S1) %*% t(S1) %*% mutable_weight %*% t(mutable_base)
  new_index <- c(determined, new_basis)
  reconciled_y <- sMat %*% rbind(reconciled_mutable, t(base_forecasts[,immutable_basis,drop=FALSE]))
  t(reconciled_y[order(new_index),])
}
