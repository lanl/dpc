\name{plot.dpc}
\alias{plot.dpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Deep Process Convolution Posterior Means
}
\description{
This function plots the posterior means and 95\% credible intervals for the functions modeled by a deep process convolution (see Lawrence, et al.) using MCMC output obtained by using the function \code{\link{dpc}}.  This function can also return the posterior realizations of the mean used to calculate the posterior mean and credible intervals.  Alternately, it can also just return the values of the posterior mean at specific points for specific functions.
}
\usage{
\method{plot}{dpc}(output = NULL, funcnums = NULL, ngrid = 50, gridpoints = NULL,
         xlim = NULL, ylim = NULL, print_warn = TRUE, graph = TRUE,
         return_draws = FALSE, return_mean = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{An object of class \code{dpc} obtained from the \code{\link{dpc}} function.}
  \item{funcnums}{A vector of integers specifying which mean functions are to be plotted.  The mean functions are numbered in the order in which they were placed in the argument of \code{Y} in the \code{\link{dpc}} function.  By default, all the functions are plotted.  The value of \code{funcnums} is sorted in ascending order in the function.}
  \item{ngrid}{The number of gridpoints to plot the posterior mean.  Larger numbers will make the plot look smoother, but will take longer.  Defaults to 50.}
  \item{gridpoints}{Numeric vector specifying which locations the posterior mean should be plotted/evaluated.}
  \item{xlim}{A vector of length 2 indicating the limits of the x-axis. Default is determined from the range of the locations supplied in the original data to \code{\link{dpc}}.}
  \item{ylim}{A vector of length 2 indicating the limits of the y-axis. Default is determined from the realized values of the mean function.}
  \item{print_warn}{Logical.  If \code{TRUE}, warning messages in the C code are displayed.}
  \item{graph}{Logical.  If \code{TRUE}, the mean functions specified in \code{funcnums} are plotted, otherwise no plotting is performed.}
  \item{return_draws}{Logical.  If \code{TRUE}, the MCMC draws of the posterior mean are returned in an array with three dimensions indexing the gridpoints, MCMC iterations, and the ordered function numbers of \code{funcnums}, respectively.}
  \item{return_mean}{Logical.  If \code{TRUE}, the mean functions of the functions in \code{funcnums} are returned at the values specified in \code{gridpoints}.}
}
\value{By default, nothing is returned. If \code{return_draws=TRUE}, a list of length 2 is returned containing the grid locations (\code{grid}) and MCMC draws of the posterior mean in an array with three dimensions indexing the gridpoints, MCMC iterations, and the ordered function numbers of \code{funcnums}, respectively.  If \code{return_mean=TRUE}, then a list of length 2 is returned with the grid locations (\code{grid}) and a matrix (\code{pmean}) with each row corresponding to the mean function at the specified gridpoints; the first column gives the function number.  If both are \code{TRUE}, a list containing the grid, posterior draws and mean functions are returned labeled \code{grid}, \code{pdraws} and \code{pmean}, respectively.  
}
\references{
Lawrence, Earl, et al. "The coyote universe. III. simulation suite and precision emulator for the nonlinear matter power spectrum." The Astrophysical Journal 713.2 (2010): 1322.
}
\author{
Richard D. Payne
}
\seealso{
\code{\link{dpc}}, \code{\link{getOutput}}
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

# Plot all three posterior mean functions and credible intervals
par(mfrow=c(3,1))
plot(out)

# Plot just one of the functions
par(mfrow=c(1,1))
plot(out,funcnums=1)

# Plot first and third posterior means and credible intervals
par(mfrow=c(2,1))
plot(out,funcnums=c(1,3))
}
}
