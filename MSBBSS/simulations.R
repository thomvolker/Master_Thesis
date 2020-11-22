
## Load required packages
library(MASS) # for drawing multivariate normal data
library(tidyverse)
library(magrittr)
library(furrr)

## Set seed 
seed <- as.integer(123)

make_coefs <- function(r2, ratio, rho) {
  
  # specify the predictor matrix, including the relative contribution of
  # all combinations
  beta_mat <- ratio %*% t(ratio)
  beta_mat[upper.tri(beta_mat, diag = TRUE)] <- 0
  # factor out the common term in all coefficients
  b <- sqrt((r2 / (sum(ratio^2) + 2 * sum(beta_mat * rho))))
  beta <- b * ratio
  beta
}

## sample data
gen_norm <- function(r2, betas, rho, n) {
  X <- mvrnorm(n = n, mu = rep(0, length(betas)), Sigma = rho)
  y <- X %*% betas + rnorm(n, 0, sqrt(1 - r2))
  return(bind_cols(Y = y, as.data.frame(X)))
}

## Specify relative importance of the regression coefficients
ratio_beta <- c(1,2,3,4)
## Specificy covariance matrix of the regression coefficients 
## (equals the correlation matrix of the regression coefficients)
rho <- diag(length(ratio_beta))
rho[rho != 1] <- 0.2

## r2 of the regression model
r2 <- c(.02, .09, .25)

## specify the sample size
n <- seq(25, 500, by = 25)

## number of simulations 
nsim <- 100

plan(multisession)


## Now, create a matrix containing the possible combinations of 
## the conditions
conditions <- expand_grid(n = n, r2 = r2) %>%
  mutate(betas = map(r2, ~ make_coefs(., ratio_beta, rho)))

data <- conditions %>%
  mutate(data  = pmap(., function(r2, betas, n) gen_norm(r2, betas, rho, n)))
  
data %>%
  mutate(model = map(data, ~ lm(Y ~ ., .))) %>%
  mutate(R2 = map(model, broom::glance)) %>%
  unnest(R2) %>%
  group_by(r2) %>%
  summarize(r.squared = mean(r.squared))



