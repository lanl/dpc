.get_TALLDAT_As <-
function(Y,vars){
  # Get size of data
  N <- length(Y)
  NIs <- numeric(N)
  NN <- 0 # total number of observations
  for(i in 1:N){
    ni <- length(Y[[i]])
    NIs[i] <- ni
    for(j in 1:ni){
      NN <- NN + nrow(Y[[i]][[j]])
    }
  }
  TALLDAT <- matrix(NA,NN,5)
  
  
  # get all unique locations
  locs <- NULL
  # TALLDAT <- NULL # Needed to pass easily to C
  startind <- 1
  for(i in 1:N){
    # number of observed reps for each smooth function i
    # print(paste0(i,", ",ni))
    # thetimes <- numeric(6)
    for(j in 1:NIs[i]){
      # end0 <- proc.time()
      dat <- as.matrix(Y[[i]][[j]])
      # end1 <- proc.time()
      # Sort data
      tempind <- order(dat[,1])
      dat <- dat[tempind,]
      # end2 <- proc.time()
      Y[[i]][[j]] <- dat
      # end3 <- proc.time()
      # ASSUME that location is the first column
      locs <- unique(c(dat[,1],locs))
      # end4 <- proc.time()
      dat <- cbind(dat,1/vars[[i]][[j]][tempind],i-1,j-1)  # Correspond to a C index
      if(length(vars[[i]][[j]]) != length(tempind)){
        stop("Variances should have same length as data")
      }
      # end5 <- proc.time()
      # TALLDAT <- rbind(TALLDAT,dat)
      n_ij <- nrow(dat)
      TALLDAT[startind:(startind + n_ij - 1),] <- dat
      startind <- startind + n_ij
      # end6 <- proc.time()
      
      # for(k in 1:6){
      #   thetimes[k] <- thetimes[k] + (get(paste0("end",k)) -
      #                                   get(paste0("end",k-1)))[3]
      # }
      
    }
    # print(sprintf("%.2f %.2f %.2f %.2f %.2f %.2f",thetimes[1],
    #        thetimes[2],thetimes[3],thetimes[4],thetimes[5],thetimes[6]))
  }
  # Find all unique locations for P_i
  locs <- sort(unique(locs))
  # Create the projection matrices, A_ij (simplifies to an index)
  As <- list(list())
  Avec <- numeric(nrow(TALLDAT))
  startind <- 1
  for(i in 1:N){
    Astmp <- list()
    # thetimes <- numeric(5)
    for(j in 1:NIs[i]){
      # end0 <- proc.time()
      Ysub <- Y[[i]][[j]][,1]
      # end1 <- proc.time()
      atmp <- numeric(length(Ysub))
      # end2 <- proc.time()
      for(k in 1:length(Ysub)){
        atmp[k] <- which(Ysub[k] == locs) - 1 # -1 for C index
      }
      # end3 <- proc.time()
      # Astmp[[j]] <- which(locs %in% Y[[i]][[j]][,1])
      Astmp[[j]] <- atmp
      # end4 <- proc.time()
      n_ij <- length(Ysub)
      Avec[startind:(startind + n_ij -1)] <- Astmp[[j]]
      startind <- startind + n_ij
      # Avec[TALLDAT[,4] == (i-1) & TALLDAT[,5] == (j-1)] <- Astmp[[j]]
      
      # end5 <- proc.time()
      # for(k in 1:5){
      #   thetimes[k] <- thetimes[k] + (get(paste0("end",k)) -
      #                                  get(paste0("end",k-1)))[3]
      # }
    }
    # print(sprintf("%.2f %.2f %.2f %.2f %.2f",thetimes[1],
    #         thetimes[2],thetimes[3],thetimes[4],thetimes[5]))
    As[[i]] <- Astmp
  }
  TALLDAT <- cbind(TALLDAT,Avec)

  output <- list(TALLDAT=TALLDAT,As=As,locs=locs,N=N,NIs=NIs)
}
