.udraw <-
function(Y,vars,poutput,outdir=NULL,
                 alpha_u=1,alpha_v=1,
                 beta_u=.001,beta_v=.001,
                 seed=NULL,print_warn=1,returnUs=FALSE){
  message("NOTE: Drawing new Us.  May replace old draws (if they exist).")
  if(is.null(outdir) & !returnUs){
    stop("Must specify 'outdir' or returnUs=TRUE")
  }
  
  res <- .get_TALLDAT_As(Y,vars)
  N <- res$N
  NIs <- res$NIs
  TALLDAT <- res$TALLDAT
  As <- res$As
  locs <- res$locs
  rm(res); gc()
  
  V <- poutput$vs
  delta <- poutput$postdraws$delta
  tau2_u <- poutput$postdraws$tau2_u
  tau2_v <- poutput$postdraws$tau2_v
  k_u <- poutput$k_u
  k_v <- poutput$k_v
  
  
  n_klocs <- length(locs)
  n_u <- length(k_u)
  n_v <- length(k_v)
  ni_max <- max(NIs)
  nmcmc <- length(delta)
  ntalldat <- nrow(TALLDAT)
  
  if(is.null(seed)){
    seed <- unclass(Sys.time()) # Obtain random seed
  }

  UU <- matrix(0,nmcmc*N,n_u + 2)
  
  # dyn.load("udraw.so")
  status <- 0
  out <- .C("udraw",as.double(alpha_u),as.double(alpha_v),as.integer(TALLDAT[,6]),
            as.double(beta_u),as.double(beta_v),
            as.double(delta),as.integer(TALLDAT[,4]),as.integer(TALLDAT[,5]),as.double(TALLDAT[,1]),
            as.double(k_u),as.double(k_v),as.double(locs),
            as.integer(N),as.integer(n_klocs),as.integer(n_u),as.integer(n_v),
            as.integer(NIs),as.integer(ni_max),as.integer(nmcmc),
            as.integer(ntalldat), as.double(TALLDAT[,2]), as.double(TALLDAT[,3]), as.integer(seed),
            as.double(tau2_u),as.double(tau2_v),as.double(UU),
            as.double(V),as.integer(print_warn),as.integer(status),PACKAGE="dpc")
  status <- out[[29]]
  if(status){
    stop("Error in 'udraw.c'")
  }
  
  Us <- matrix(out[[26]],nmcmc*N,n_u + 2,byrow=TRUE)
  Us2 <- list()
  for(i in 1:N){
    Us2[[sprintf("u%04d",i)]] <- Us[seq(from=i,to=nmcmc*N,by=N),3:(n_u + 2)]
  }
  
  if(!is.null(outdir)){
    if(substring(outdir,nchar(outdir)) != "/"){ # Add directory slash if needed
      outdir <- paste0(outdir,"/")
    }
    # write.table(Us,paste0(outdir,"Us.csv"),sep=',',row.names=FALSE,col.names=FALSE)
    udir <- paste0(outdir,"us/")
    if(!dir.exists(udir)){
      dir.create(udir,recursive=TRUE)
    }
    for(i in 1:N){
      # write.table(Us[seq(from=i,to=nmcmc*N,by=N),3:(n_u + 2)],
      #             paste0(udir,sprintf("u%04d",i,".csv")),
      #             row.names=FALSE,sep=',',col.names=FALSE)
      write.table(Us2[[i]],paste0(udir,sprintf("u%04d",i,".csv")),
                              row.names=FALSE,sep=',',col.names=FALSE)
    }
  }
  if(returnUs){
    return(Us2)
  }
}
