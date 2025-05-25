\name{get_Ksigma}
\alias{get_Ksigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct K-sigma Matrix
}
\description{
This function constructs the K-sigma matrix for deep process convolutions.
}
\usage{
get_Ksigma(locs, k_u, sigmas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{locs}{The locations at which the mean function is to be evaluated.}
  \item{k_u}{The locations of the latent variable u.}
  \item{sigmas}{The bandwidth parameter at each location, \code{locs}}
}
\details{
The matrix K-sigma is the matrix that is used to construct the mean function as K-sigma * u where u is the latent variable in the deep process convolution model.
}
\value{
A matrix of dimension \code{length(locs)} by \code{length(k_u)} where the element i,j is specified as 1/sqrt(2*pi*sigma_i^2) * exp((\code{locs[i]} - \code{k_u[j]})^2 /(2*sigma_i^2)).
}
\references{
Lawrence, Earl, et al. "The coyote universe. III. simulation suite and precision emulator for the nonlinear matter power spectrum." The Astrophysical Journal 713.2 (2010): 1322.
}
\author{
Richard D. Payne
}
\seealso{
\code{\link{get_Kdelta}}
}