# dpc: Bayesian Deep Process Convolutions

The `dpc` package provides tools to fit **Bayesian Deep Process Convolution** models for estimating nonstationary mean functions from noisy, multiresolution data. It implements the two-layer convolution approach introduced in:

> Moran, Kelly R., et al. *Bayesian “Deep” Process Convolutions: An Application in Cosmology*, 2024  
> [arXiv:2411.14747](https://arxiv.org/abs/2411.14747)

This model is particularly useful when:
- The function's smoothness varies over its domain.
- You have multiple related functions (e.g., cosmologies) sharing structural features.
- You're working with data of varying resolution and uncertainty.

## Features

- Two-layer convolution (deep process convolution)
- Learns spatially varying smoothness automatically
- Supports multiple processes sharing information
- Handles non-uniform grids and unequal variances
- Parallelized C backend for efficient MCMC

## Installation

This package depends on `RcppGSL` and requires compilation. To install from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install dpc from GitHub
devtools::install_github("lanl/dpc")
```

Ensure you have a C++ compiler and GSL (GNU Scientific Library) installed on your system. On Linux, for example:

```bash
sudo apt-get install libgsl-dev
```

On macOS:

```bash
brew install gsl
```

## Usage

### Quick Start

Here's a minimal example using a noisy version of the function `sin(1/x)`:

```r
library(dpc)

# Define function and noise
f <- function(x) sin(1/x)
sig <- function(x, c=0.5) c * rnorm(length(x), 0, sqrt(x))

# Generate toy data
y = vars = list()
y[[1]] = vars[[1]] = list()
set.seed(123)
nrep = 10; n = 100
for(jj in 1:nrep){
  x <- runif(n, 0.1, 1)
  sigx <- sig(x)
  y[[1]][[jj]] = cbind(x, f(x) + sigx)
  vars[[1]][[jj]] = sigx^2
}

# Fit the model
res = dpc(Y = y, vars = vars, nmcmc = 1000, burn = 500, num_threads = 1)

# Plot estimated mean functions with uncertainty
plot.dpc(res)
```

### Accessing Results

```r
str(res$hyperparms)  # Initial values and MCMC settings
str(res$locs)        # Unique observed locations
str(res$k_u)         # Grid for latent u process
str(res$k_v)         # Grid for latent v process
str(res$lpost)       # Posterior log-density
str(res$postdraws)   # Posterior samples of hyperparameters
str(res$us)          # Posterior draws for latent process u
str(res$vs)          # Posterior draws for latent process v
```

### More Complex Example

Estimate multiple related functions with differing noise levels:

```r
library(dpc)
set.seed(424)

# Create data from 3 functions with different smoothness and noise
n <- 75
f1 <- function(x) sin(x)
f2 <- function(x) cos(x)
f3 <- function(x) 1 + sin(x)
sds <- c(0.5, 1, 1.5)

Y <- vars <- list()
for(i in 1:3){
  Ytmp <- varstmp <- list()
  for(j in 1:5){
    x <- sort(runif(n, 0, 2*pi))
    y <- do.call(paste0("f", i), list(x = x)) + rnorm(n, 0, sds[i])
    Ytmp[[j]] <- cbind(x, y)
    varstmp[[j]] <- rep(sds[i]^2, n)
  }
  Y[[i]] <- Ytmp
  vars[[i]] <- varstmp
}

# Fit model
out <- dpc(Y, vars, nmcmc = 1000, burn = 1000, n_u = 15, n_v = 5, seed = 1, num_threads = 1, draw_u = TRUE)

# Plot results
plot.dpc(out)
```

## Reference

Please cite the following if you use `dpc` in your work:

```
@article{moran2024bayesian,
  title={Bayesian “Deep” Process Convolutions: An Application in Cosmology},
  author={Moran, Kelly R and Payne, Richard D and others},
  journal={arXiv preprint arXiv:2411.14747},
  year={2024}
}
```

## License

This package was approved for Open-Source Assertion as O4860. Please see LICENSE for further copyright details.

## Contact

For issues, questions, or contributions, please contact: Kelly R. Moran at krmoran at lanl dot gov.
