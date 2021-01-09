
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

gen_log <- function(r2, betas, rho, n) {
  dat <- 
    gen_norm(r2, betas, rho, n) %>%
    mutate(Y = rbinom(n, 1, exp(Y) / (1 + exp(Y))))
  
  dat
}

gen_prob <- function(r2, betas, rho, n) {
  dat <- 
    gen_norm(r2, betas, rho, n) %>%
    mutate(Y = rbinom(n, 1, pnorm(Y)))
  
  dat
}

get_bf <- function(bf) {
  if (class(bf)[2] != "bain") stop("Input must have class bain")
  bf$fit$BF.u
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

n <- doubling_seq(start = 25, n = 6)

## number of simulations 
nsim <- 500

## Now, create a matrix containing the possible combinations of 
## the conditions
conditions <- expand_grid(nsim = 1:nsim, n = n, r2 = r2) %>%
  mutate(betas = map(r2, ~ make_coefs(., ratio_beta, rho)))

hypo <- paste0("V1<V2<V3<V4;V1=V2=V3=V4")

plan(multisession)

datasets <- conditions %>%
  mutate(lin_dat  = future_pmap(., function(nsim, r2, betas, n) gen_norm(r2, betas, rho, n), .progress = TRUE, .options = future_options(seed = as.integer(123))),
         log_dat  = future_pmap(., function(nsim, r2, betas, n) gen_log(r2, betas, rho, n), .progress = TRUE, .options = future_options(seed = as.integer(124))),
         prob_dat = future_pmap(., function(nsim, r2, betas, n) gen_prob(r2, betas, rho, n), .progress = TRUE, .options = future_options(seed = as.integer(125))))

models <- datasets %>%
  mutate(lin_mod  = map(lin_dat,    ~lm(Y ~ ., data = .), .progress = TRUE),
         log_mod  = map(log_dat,   ~glm(Y ~ ., data = ., family = binomial(link = "logit"))),
         prob_mod = map(prob_dat,  ~glm(Y ~ ., data = ., family = binomial(link = "probit"))))

bain_out <- models %>%
  mutate(bf_lin   = map(lin_mod,  ~bain(., hypothesis = hypo)),
         bf_log   = map(log_mod,  ~bain(., hypothesis = hypo)),
         bf_prob  = map(prob_mod, ~bain(., hypothesis = hypo)))




test_stats <- bain_out %>%
  mutate(lin_t = map(lin_mod, function(x) x %>% summary %>% coef %>% .[,3]),
         log_z = map(log_mod, function(x) x %>% summary %>% coef %>% .[,3]),
         prob_z = map(prob_mod, function(x) x %>% summary %>% coef %>% .[,3])) %>%
  unnest(c(lin_t, log_z, prob_z)) %>%
  select(nsim, n, r2, lin_t, log_z, prob_z) %>%
  unnest(c(lin_t, log_z, prob_z), .id = "ID") %>%
  mutate("var" = rep(c("Int", "V1", "V2", "V3", "V4"), nrow(conditions))) %>%
  group_by(n, r2, var) %>%
  summarize(lin_t = mean(lin_t),
            log_z = mean(log_z),
            prob_z = mean(prob_z))

test_stats %>% 
  pivot_longer(cols = c(lin_t, log_z, prob_z), names_to = "model") %>%
  ggplot(aes(x = n)) +
  geom_line(aes(y = value, col = var)) +
  facet_wrap(r2 ~ model) +
  theme_minimal()


bain_out %>%
  mutate(bf_lin_h1  = map_dbl(bf_lin, function(x) get_bf(x)[1]),
         bf_log_h1  = map_dbl(bf_log, function(x) get_bf(x)[1]),
         bf_prob_h1 = map_dbl(bf_prob, function(x) get_bf(x)[1]),
         bf_h1      = bf_lin_h1 * bf_log_h1 * bf_prob_h1) %>%
  group_by(n, r2) %>%
  summarize(bf_lin_h1 = mean(bf_lin_h1),
            bf_log_h1 = mean(bf_log_h1),
            bf_prob_h1 = mean(bf_prob_h1),
            BF_H1 = mean(bf_h1)) %>%
  pivot_longer(cols = c(bf_lin_h1, bf_log_h1, bf_prob_h1)) %>%
  ggplot(aes(x = n, y = value, col = as.factor(r2))) +
  facet_wrap(~ name) +
  geom_line() +
  geom_abline(intercept = 1, slope = 0) +
  theme_minimal()

output$lin_mod[[3]] %>% summary
output <- output %>%
  mutate(pmp = map(bfs, function(x) x$fit$PMPb[1]),
         pmpe = map(bfs, function(x) x$fit$PMPb[2])) %>%
  unnest(c(pmp, pmpe))

conditions[900,]
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

## to do next
## create function that can combine multiple bain objects and returns PMPs
## and bayes factors

