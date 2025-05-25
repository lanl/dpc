get_Kdelta <-
function(locs,k_v,delta){
  nloc <- length(locs)
  n_v <- length(k_v)
  K_delta <- matrix(NA,nloc,n_v)
  for(j in 1:n_v){
    K_delta[,j] <- dnorm(locs,k_v[j],delta)
  }
  K_delta <- K_delta / rowSums(K_delta)
  return(K_delta)
}
