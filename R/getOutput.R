getOutput <-
function(path="./",draw_u=FALSE,
         Y=NULL,vars=NULL,seed=NULL,print_warn=TRUE,lite=FALSE){
 
  if(draw_u){
    if(is.null(Y)) stop("Must provide 'Y' argument when drawing Us")
    if(is.null(vars)) stop("Must provide 'vars' argument when drawing Us")
  }
  
  if(substr(path,nchar(path),nchar(path)) != '/'){
    path = paste0(path,"/")
  }
  
  k_u <- scan(paste0(path,"k_u.txt"),quiet=TRUE)
  k_v <- scan(paste0(path,"k_v.txt"),quiet=TRUE)
  locs <- scan(paste0(path,"k_locs.txt"),quiet=TRUE)
  lpost <- scan(paste0(path,"log_post.txt"),quiet=TRUE)
  out <- read.csv(paste0(path,"output.csv"))
  vs <- as.matrix(read.csv(paste0(path,"v_output.csv"),header=FALSE))
  hyperparms <- read.csv(paste0(path,"hyperparms.csv"))
  OUT <- list(locs=locs,k_u=k_u,k_v=k_v,lpost=lpost,
              postdraws=out,vs=vs,hyperparms=hyperparms)
  
  us <- NULL
  # fname_u <- paste0(path,"Us.csv")
  udir <- paste0(path,'us/')
  if(dir.exists(udir)){
    if(!lite){ # read in all u files
      us <- list()
      ufiles <- dir(udir)
      for(i in 1:length(ufiles)){
        us[[paste0("u",sprintf("%04i",i))]] <- as.matrix(read.csv(paste0(udir,ufiles[i]),header=FALSE))
      }
    }else{
      us <- udir
    }
    OUT[["us"]] <- us
  }else if(draw_u){
    us <- .udraw(Y,vars,OUT,outdir=path,
          alpha_u=OUT$hyperparms$alpha_u,
          alpha_v=OUT$hyperparms$alpha_v,
          beta_u=OUT$hyperparms$beta_u,
          beta_v=OUT$hyperparms$beta_v,
          seed=seed,print_warn=print_warn,
          returnUs=TRUE)
    OUT[["us"]] <- us
  }else if(!dir.exists(udir)){
    message("NOTE: Directory for 'us' cannot be found. Consider drawing them.")
  }
  OUT[["lite"]] <- lite
  class(OUT) <- "dpc"
  return(OUT)
}
