dpc <- function(Y,vars,nmcmc=10000,burn=1000,thin=1,n_u=100,n_v=10,
                       k_u=NULL,k_v=NULL,
                       v=NULL,tau2_u=.01,tau2_v=1,delta=.5,delta_max=NULL,
                       alpha_u=1,alpha_v=1,beta_u=.001,beta_v=.001,
                       sprop_v=NULL,sprop_tau2_u=.05,sprop_delta=.05,
                       n_print=100,out_dir="./output/",seed=0,draw_u=TRUE,
                       num_threads=0,print_warn=1,lite=FALSE){
  print("Preparing Data...",quote=FALSE)
  
  res <- .get_TALLDAT_As(Y,vars)
  N <- res$N
  NIs <- res$NIs
  TALLDAT <- res$TALLDAT
  As <- res$As
  locs <- res$locs
  rm(res); gc()
  
  max_ni <- max(NIs)
  ntalldat <- nrow(TALLDAT)
  n_klocs <- length(locs)
  # Validate ulim and vlim values
  # if(is.null(ulim) || length(ulim) != 2){
  #   if(!is.null(ulim)){
  #     stop("ulim must be a vector of length 2 or NULL")
  #   }
  ulim=c(min(locs),max(locs))
  #}
  # if(is.null(vlim) || length(vlim) != 2){
  #   if(!is.null(vlim)){
  #     stop("vlim must be a vector of length 2 or NULL")
  #   }
  vlim <- c(min(locs),max(locs))
  #}
  
  if(is.null(k_u)){
    bnds <- .gridminmax(ulim[1],ulim[2],n_u)
    k_u <- seq(from=bnds[1],to=bnds[2],length.out=n_u)
  }else{
    n_u <- length(k_u)
    message(paste0("NOTE: n_u overidden and set to n_u=",n_u))
  }
  if(is.null(k_v)){
    bnds <- .gridminmax(vlim[1],vlim[2],n_v)
    k_v <- seq(from=bnds[1],to=bnds[2],length.out=n_v)
  }else{
    n_v <- length(k_v)
    message(paste0("NOTE: n_v overidden and set to n_v=",n_v))
  }
  if(is.null(delta_max)){
    delta_max <- vlim[2] - vlim[1] # Ensure delta_max isn't too small
  }
  if(delta > delta_max) stop("Delta must be less than delta_max")
  
  
  # Initial values
  if(is.null(v)){
    v <- rep(.1,n_v)
  }else if(length(v) != n_v){
    stop("v must have length n_v.")
  }
  
  # Proposal distributions
  if(is.null(sprop_v)){
    sprop_v <- rep(.1,n_v)
  }else if(length(sprop_v) != n_v){
    stop("sprop_v must have length n_v.")
  }
  # sprop_tau2_u <- .2
  # sprop_delta <- .05
  # n_print <- 100
  # 
  # add slash at end if needed
  if(substr(out_dir,nchar(out_dir),nchar(out_dir)) != '/'){
    out_dir = paste0(out_dir,"/")
  }
  # Create directory if it doesn't exist
  if(!dir.exists(out_dir)){
    dir.create(out_dir,recursive=TRUE)
  }
  if(draw_u){
    udirec <- paste0(out_dir,"us/")
    if(!dir.exists(udirec)){
      dir.create(udirec,recursive=TRUE)
    }
  }
  
  # Save starting values
  start_v <- v
  start_tau2_u <- tau2_u
  start_tau2_v <- tau2_v
  start_delta <- delta
  
  dpc_status <- 0
  
  # Call C function
  nothing <- .C('dpc',as.double(TALLDAT[,1]),as.double(TALLDAT[,2]),
                as.double(TALLDAT[,3]),
     as.integer(TALLDAT[,4]),as.integer(TALLDAT[,5]),
     as.integer(TALLDAT[,6]),as.double(locs),as.integer(NIs),
     as.integer(N),as.integer(ntalldat),as.double(v),as.integer(n_klocs),
     as.integer(n_u),as.integer(n_v),as.double(k_u),as.double(k_v),
     as.double(tau2_u),as.double(tau2_v),as.double(delta),
     as.double(alpha_u),as.double(alpha_v),as.double(beta_u),
     as.double(beta_v),as.integer(nmcmc),as.integer(burn),
     as.double(sprop_v),as.double(sprop_tau2_u),
     as.double(sprop_delta),as.integer(n_print),
     as.character(out_dir),as.integer(seed),as.integer(draw_u),
     as.integer(num_threads),as.double(delta_max),
     as.integer(dpc_status),as.integer(print_warn),
     as.integer(thin),PACKAGE="dpc")
  
  dpc_status <- nothing[[35]]
  if(dpc_status > 0){
    return(NULL);
  }
  
  # Print Diagnostics
  cat(readLines(paste0(out_dir,"mcmc_info.txt")),sep="\n")
  
  # Save hyperparameters to file
  hyperparms <- data.frame(alpha_u=alpha_u,alpha_v=alpha_v,beta_u=beta_u,beta_v=beta_v,delta_max=delta_max,N=N)
  write.table(hyperparms,paste0(out_dir,"hyperparms.csv"),sep=',',row.names=FALSE)

  # Read back in data for user to use
  out <- getOutput(out_dir,lite=lite)
  
  
  return(out)
  
  # lpost <- scan(paste0(out_dir,"log_post.txt"),quiet=TRUE)
  # out <- read.csv(paste0(out_dir,"output.csv"))
  # vs <- read.csv(paste0(out_dir,"v_output.csv"),header=FALSE)
  # us <- NULL
  # if(draw_u){
  #   us <- read.csv(paste0(out_dir,"Us.csv"),header=FALSE)
  # }
  # return(
  #   list(locs=locs,k_u=k_u,k_v=k_v,lpost=lpost,
  #        postdraws=out,vs=vs,us=us)
  # )
}









