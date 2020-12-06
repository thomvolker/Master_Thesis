
## Load required packages
library(MASS) # for drawing multivariate normal data
library(tidyverse)
library(magrittr)
library(furrr)
library(bain)

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
doubling_seq <- function(start, end = NULL, n = NULL) {
  if (is.null(end) && is.null(n)) stop("Either `end` or `n` must be specified")
  if (is.null(end)) return(start * 2^(0:(n-1)))
  if (is.null(n)) return(start * 2^(0:log2(end/start)))
}

n <- doubling_seq(start = 25, n = 15)

## number of simulations 
nsim <- 50

plan(sequential)


## Now, create a matrix containing the possible combinations of 
## the conditions
conditions <- expand_grid(nsim = 1:nsim, n = n, r2 = r2) %>%
  mutate(betas = map(r2, ~ make_coefs(., ratio_beta, rho)))

hypo <- paste0("V1<V2<V3<V4;V1=V2=V3=V4")


output <- conditions %>%
  mutate(data  = pmap(., function(nsim, r2, betas, n) gen_norm(r2, betas, rho, n)),
         model = map(data,  ~lm(Y ~ ., .)),
         bfs   = map(model, ~bain(., hypothesis = hypo)))

output <- output %>%
  mutate(pmp = map(bfs, function(x) x$fit$PMPb[1]),
         pmpe = map(bfs, function(x) x$fit$PMPb[2])) %>%
  unnest(c(pmp, pmpe))


output %>%
  group_by(n, r2) %>%
  summarize(pmp = mean(pmp),
            pmpe = mean(pmpe)) %>%
  ggplot(mapping = aes(x = n, 
                     y = pmp, 
                     col = as.factor(r2), 
                     group = as.factor(r2))) +
    geom_line() +
    geom_line(aes(y = pmpe), linetype = "dashed") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")


data %>%
  rowid_to_column(var = "ID") %>%
  mutate(model = map(data, ~ lm(Y ~ ., .))) %>%
  mutate(stats = map(model, broom::glance)) %>%
  unnest(stats) %>%
  group_by(r2) %>%
  summarize(r2 = mean(adj.r.squared))

         