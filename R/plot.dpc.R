plot.dpc <- function (output = NULL, funcnums = NULL, ngrid = 50, 
                      gridpoints = NULL,
                      xlim = NULL, ylim = NULL, print_warn = TRUE, graph = TRUE,
                      return_draws = FALSE,return_mean = FALSE) {
  
  if(is.null(output)) stop("Must provide 'output' argument.")
  # if((return_draws | return_mean) & is.null(gridpoints)){
  #   stop("Must specify gridpoints if return_draws or return_mean is TRUE.")
  # }
  
  if(is.null(output$us)) stop("Must have u in order to plot mean function.")
  
  
  # Add slash to end of path if needed
  # if(!is.null(path)){
  #   if(substr(path,nchar(path),nchar(path)) != '/'){
  #     path = paste0(path,"/")
  #   }
  # }
  # 
  # if (!is.null(output) && is.null(path)) {
  # }else if (is.null(output) && !is.null(path)) {
  #   output <- getOutput(path)
  # }else if(!is.null(output) && !is.null(path)){ # must be lite
  #    if(!output$lite){
  #      stop("You may only specify 'output' and 'path' if output is 'lite'.")
  #    }
  # }else{
  #   error("Must provide either the output or the path.")
  # }
  
  # funcnums_all <- unique(output$us[, 2])
  # N <- length(funcnums_all)
  N <- output$hyperparms$N
  funcnums_all <- 1:N
  nmcmc <- length(output$lpost)
  n_u <- length(output$k_u)
  if (is.null(funcnums)) {
    funcnums <- funcnums_all
  }else if (!all(funcnums %in% funcnums_all)) {
    stop(paste0("funcnums must be between 1 and ", max(funcnums_all)))
  }
  funcnums <- sort(funcnums)
  if (is.null(xlim)) {
    xlim <- c(min(output$locs), max(output$locs))
  }
  ylim_orig <- ylim
  if(is.null(gridpoints)){
    thegrid <- seq(from = min(output$locs), to = max(output$locs), 
                   length.out = ngrid)
  }else{
    thegrid <- sort(gridpoints)
    ngrid <- length(gridpoints)
  }
  pdraws <- array(0,dim=c(ngrid,nmcmc,length(funcnums)))
  # for (i in 1:nmcmc) {
  #   K_delta <- get_Kdelta(thegrid, output$k_v, output$postdraws$delta[i])
  #   sigmas <- c(K_delta %*% output$vs[i, ])
  #   K_sigma <- get_Ksigma(thegrid, output$k_u, sigmas)
  #   for (j in 1:length(funcnums)) {
  #     funcnum <- funcnums[j]
  #     if(all(sigmas > 0)){
  #       pdrawval <- K_sigma %*% output$us[funcnum + (i - 1) * N, 3:(n_u + 2)]
  #     }else{
  #       pdrawval <- NA
  #     }
  #     pdraws[, i, j] <- pdrawval
  #   }
  # }
  n_func <- length(funcnums)
  n_klocs <- length(output$locs)
  n_u <- length(output$k_u)
  n_v <- length(output$k_v)
  status <- 0
  
  if(output$lite){  # Read in only necessary Us
    udir <- output$us
    funcnums_lite <- 1:length(funcnums)
    US <- matrix(NA,nmcmc*length(funcnums),n_u)
    for(i in 1:length(funcnums)){
      US[seq(from=i,to=nmcmc*length(funcnums),by=length(funcnums)),] <- as.matrix(read.csv(paste0(udir,sprintf("u%04d.csv",funcnums[i])),header=FALSE))
    }
    funcnums_c <- funcnums_lite
    N_c <- length(funcnums)
  }else{
    US <- matrix(NA,nmcmc*N,n_u)
    for(i in 1:N){
      US[seq(from=i,to=nmcmc*N,by=N),] <- output$us[[sprintf("u%04d",i)]]
    }
    funcnums_c <- funcnums
    N_c <- N
  }
  
  res <- .C("postdraws",as.double(pdraws),as.double(output$postdraws$delta),
                as.integer(funcnums_c),
     as.double(output$k_u),as.double(output$k_v),as.double(output$locs),
     as.integer(n_func),as.integer(N_c),as.integer(ngrid),as.double(thegrid),
     as.integer(nmcmc),
     as.integer(n_klocs),as.integer(n_u),as.integer(n_v),
     as.double(US),as.double(output$vs),as.integer(print_warn),
     as.integer(status),
     PACKAGE="dpc")
  status <- res[[18]]
  if(status){
    stop("Error occurred in C function 'postdraws()'")
  }
  pdraws <- array(res[[1]],dim=c(ngrid,nmcmc,n_func))
  ind <- which(pdraws == -999,arr.ind = TRUE)
  pdraws[ind] <- NA
  
  PMEAN <- NULL
  if(graph | return_mean){
    PMEAN <- matrix(NA,length(funcnums),ngrid)
    for(j in 1:length(funcnums)){
      PMEAN[j,] <- apply(pdraws[,,j], 1, mean,na.rm=TRUE)
    }
  }

  if(graph){
    for(j in 1:length(funcnums)){
      if (is.null(ylim_orig)) {
        ylim <- c(min(pdraws[,,j],na.rm=TRUE), max(pdraws[,,j],na.rm=TRUE))
      }
      # pmean <- apply(pdraws[,,j], 1, mean,na.rm=TRUE)
      pmean <- PMEAN[j,]
      plot(thegrid, pmean, col = "cornflowerblue", type = "l", 
           lwd = 1.5, xlim = xlim, ylim = ylim, 
           main = paste0("Function ",funcnums[j]),
           xlab = "", ylab = "")
      qtls <- apply(pdraws[,,j], 1, quantile, c(0.025, 0.975),na.rm=TRUE)
      lines(thegrid, qtls[1, ], lty = 2, col = "gray50")
      lines(thegrid, qtls[2, ], lty = 2, col = "gray50")
    }
  }
  if(return_draws & return_mean) return(list(pdraws=pdraws,pmean=cbind(funcnums,PMEAN),grid=thegrid))
  if(return_draws) return(list(pdraws=pdraws,grid=thegrid))
  if(return_mean){
    return_mean_dat <- cbind(funcnums,PMEAN)
    colnames(return_mean_dat)[1] <- 'func_num'
    return(list(pmean=return_mean_dat,grid=thegrid))
  } 
}
