\name{dpc}
\alias{dpc}
\title{Deep Process Convolutions}
\description{
Performs a Bayesian analysis of data using deep process convolutions (dpc)
to model the mean.
}
\usage{
dpc(Y, vars, nmcmc = 10000, burn = 1000, thin = 1, n_u = 100, n_v = 10, 
    k_u = NULL, k_v = NULL, v = NULL, tau2_u = 0.01, tau2_v = 1, 
    delta = 0.5, delta_max = NULL, alpha_u = 1, alpha_v = 1, 
    beta_u = 0.001, beta_v = 0.001, sprop_v = NULL, 
    sprop_tau2_u = 0.05, sprop_delta = 0.05, n_print = 100,
    out_dir = "./output/", seed = 0, draw_u = TRUE,
    num_threads = 0, print_warn = 1, lite = FALSE)
}
\arguments{
  \item{Y}{A list of n lists.  The ith list should consist of n_i matrices/dataframes of realizations which are assumed to have the same mean function.  Specifically, the i,jth element of \code{Y} (i.e. \code{Y[[i]][[j]]}) should be a numeric matrix or dataframe consisting of two columns.  The first column consists of the locations of the observed realization and the second column are the values of the realization (data).}
  \item{vars}{A list of n lists which contain the variances corresponding to the locations in each realization.  That is, \code{vars[[i]][[j]]} should contain a vector of the same length as \code{nrow(Y[[i]][[j]])} specifying the variance at each location.}
  \item{nmcmc}{A positive integer indicating the number of MCMC iterations to collect. Defaults to 10^4.}
  \item{burn}{A positive integer indicating the number of MCMC iterations to perform as burn in before collecting \code{nmcmc} iterations. Defaults to 10^3.}
  \item{thin}{A positive integer. After burn in, only every \code{thin} MCMC sample will be kept.  Default is 1.}
  \item{n_u}{The number of grid locations for u.  Will be overridden to be \code{length(k_u)} if \code{k_u} is not \code{NULL}. Defaults to 100.}
  \item{n_v}{The number of grid locations for v.  Will be overridden to be \code{length(k_v)} if \code{k_v} is not \code{NULL}.  Defaults to 10}
  \item{k_u}{The grid locations (numeric) for the latent variable u.  If \code{NULL}, the grid points are automatically determined using the specified locations in \code{Y}.}
  \item{k_v}{The grid locations (numeric) for the latent variable v.  If \code{NULL}, the grid points are automatically determined using the specified locations in \code{Y}.}
  \item{v}{A numeric vector of length \code{n_v} specifying the starting values for the latent function v.}
  \item{tau2_u}{A positive value for the starting value of tau^2_u. Default is .01.}
  \item{tau2_v}{A positive value for the starting value of tau^2_v. Default is 1.}
  \item{delta}{A positive value for the starting value of delta. Default is .5}
  \item{delta_max}{The maximum value of delta allowed.  Default is the range of the realized locations.}
  \item{alpha_u}{A positive value for the hyperparameter alpha_u.  Default is 1.}
  \item{alpha_v}{A positive value for the hyperparameter alpha_v.  Default is 1.}
  \item{beta_u}{A positive value for the hyperparameter beta_u. Default is .001.}
  \item{beta_v}{A positive value for the hyperparameter beta_v  Default is .001.}
  \item{sprop_v}{A numeric vector of length \code{n_v} of proposal standard deviations for the MCMC algorithm.  Default is \code{c(.1,...,.1)}.}
  \item{sprop_tau2_u}{A positive value indicating the proposal standard deviation for tau^2_u.  Default is .05.}
  \item{sprop_delta}{A positive value indicating the proposal standard deviation for delta. Default is .05.}
  \item{n_print}{An integer.  Diagnostic information of the MCMC chain will be printed every \code{n_print} iterations.}
  \item{out_dir}{The directory where the results are to be saved.  Default is "./output/"}
  \item{seed}{Sets the seed so the results are reproducible.  If 0, a random seed is chosen based on the time.  Default is 0.}
  \item{draw_u}{Logical (default \code{TRUE}) indicating if the latent variable u should be drawn (\code{TRUE}) or not (\code{FALSE}).}
  \item{num_threads}{An integer indicating how many threads should be used to compute the likelihood in parallel if OpenMP is available.  If 0, all available threads are used. Defaults to 0.}
  \item{print_warn}{An integer indicating if warnings from the C code should be printed (1) or not (0).  Defaults to 1.}
  \item{lite}{Logical indicating if the draws of the latent variable u should be returned if they are drawn.  When \code{lite=TRUE}, the values of u are not returned, but will be saved to file if \code{draw_u=1}.  This is intended to save memory for large datasets and/or long MCMC runs.}
}
\details{
The deep process convolution model is the same as that in Lawrence, with slightly different notation and generalization.  The data are not restricted to be at the same locations, but are allowed to be distinct for each realization.  For model specifics, please refer to Lawrence et al.
}
\value{
An object of class dpc.  The object is a list containing the following:

\item{hyperparms}{A dataframe with the initial hyperparmeters and number of smooth functions which were estimated.}

\item{k_u}{The grid locations of the latent vector u.}

\item{k_v}{The grid locations of the latent vector v.}

\item{lite}{A flag indicating if the draws of the latent variable u are included in the object.}

\item{locs}{A numeric vector of all of the unique locations of the data.}

\item{lpost}{A numeric vector of length \code{nmcmc} with the value of the log of the function proportional to the posterior at each iteration of the MCMC.}

\item{postdraws}{A dataframe with the posterior draws of tau^2_u, tau^2_v, and delta.}

\item{us}{If \code{lite=FALSE}, a list of draws of the latent variable u for each of the smooth functions labeled uXXXX.  If \code{lite=TRUE}, a string with the local path to the directory where the draws of u are stored.}

\item{vs}{A matrix containing the posterior draws of the latent vector v.}

}

\references{
Lawrence, Earl, et al. The coyote universe. III. simulation suite and precision emulator for the nonlinear matter power spectrum. The Astrophysical Journal 713.2 (2010): 1322.
}
\author{
Richard D. Payne
}

\seealso{
\code{\link{plot.dpc}}, \code{\link{getOutput}}
}
\examples{
\dontrun{
library(dpc)
set.seed(424)

# Create data from 3 smooth functions, 5 realizations each
vars <- Y <- list()
n <- 75 # number of observations per realization
f1 <- function(x) sin(x)
f2 <- function(x) cos(x)
f3 <- function(x) 1 + sin(x)
sds <- c(.5,1,1.5)
for(i in 1:3){
  varstmp <- Ytmp <- list()
  for(j in 1:5){
    # Random locations
    x <- sort(runif(n,0,2*pi))
    m <- do.call(paste0("f",i),alist(x=x))
    y <- m + rnorm(n,0,sds[i])
    Ytmp[[j]] <- cbind(x,y)
    varstmp[[j]] <- rep(sds[i]^2,n)
  }
  Y[[i]] <- Ytmp
  vars[[i]] <- varstmp
}

# Plot the realized data
par(mfrow=c(3,1))
for(i in 1:3){
  plot(Y[[i]][[1]],type='l')
  for(j in 2:5){
    lines(Y[[i]][[j]])
  }
}

# Run MCMC
out <- dpc(Y,vars,nmcmc=10^3,burn=10^3,n_u=15,n_v=5,seed=1,num_threads=1,draw_u = TRUE)

# Plot posterior mean functions
par(mfrow=c(3,1))
plot(out)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}% __ONLY ONE__ keyword per line
