get_Ksigma <-
function(locs,k_u,sigmas){
  nloc <- length(locs)
  n_u <- length(k_u)
  K_sigma <- matrix(NA,nloc,n_u)
  for(j in 1:n_u){
    colval <- tryCatch(dnorm(k_u[j],locs,sigmas),
                       warning=function(w) NA)
    if(length(colval) == 1){
      if(is.na(colval)) return(NA) # Invalid value
    }
    K_sigma[,j] <- colval
  }
  K_sigma <- K_sigma / rowSums(K_sigma)
  return(K_sigma)
}
