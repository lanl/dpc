.gridminmax <- function(a,b,n_u){
  n_u <- 10
  w <- 1 - 1/(2*(n_u-1))
  C <- 1 - 1/(2*(n_u-1)) - 1/(w*4*(n_u-1)^2)
  astar <- (-b/(w*2*(n_u-1)) + a) * C^(-1)
  bstar <- (-astar/(2*(n_u-1)) + b) * 1/w
  return(c(astar,bstar))
}


