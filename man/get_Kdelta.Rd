\name{get_Kdelta}
\alias{get_Kdelta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct K-delta Matrix
}
\description{
This function constructs the K-delta matrix for deep process convolutions.
}
\usage{
get_Kdelta(locs, k_v, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{locs}{The locations at which an estimate of sigma is needed.}
  \item{k_v}{The locations of the latent variable v.}
  \item{delta}{The value of delta to use in the construction of the matrix.}
}
\details{
The matrix K-delta is the matrix such that sigma = K-delta v where sigma is the bandwidth parameter for the latent variable u in the deep process convolution model.
}
\value{
A matrix of dimension \code{length(locs)} by \code{length(k_v)} where the element i,j is specified as (sqrt(2*pi*delta^2)^{-1} * exp(-(\code{locs[i]} - \code{k_v[j]})^2/(2*delta^2)).
}
\references{
Lawrence, Earl, et al. "The coyote universe. III. simulation suite and precision emulator for the nonlinear matter power spectrum." The Astrophysical Journal 713.2 (2010): 1322.
}
\author{
Richard D. Payne
}
\seealso{
\code{\link{get_Ksigma}}
}
