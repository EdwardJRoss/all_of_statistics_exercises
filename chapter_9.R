## Ex 9.2
mle_unif_mean <- function(x) {
    (min(x) + max(x)) / 2
}

unif_mle_mse <- function(a=1, b=3, n=10, trials=1000) {
    unif_samples <- matrix(runif(n*trials, a, b), nrow = trials)
    mle_estimates <- apply(unif_samples, 1, mle_unif_mean)
    mean((mle_estimates - (b+a)/2)^2)
}

unif_mle_mse(1, 3, 10, 1000)

## MSE approx proportional to 1/n^2 for large n, some fixed a, b
unif_series_size <- vector("double", 100)
for (i in seq_along(unif_series_size)) {
    unif_series_size[i] <- unif_mle_mse(0, 1, i)
}
plot(1:100, unif_series_size^(-0.5), 'l')

## MSE proportional to (b-a)^2 for large n, fixed (a+b)/2, n
unif_series_scale <- vector("double", 100)
for (i in seq_along(unif_series_scale)) {
    unif_series_scale[i] <- unif_mle_mse(0, i, 10)
}
plot(1:100, unif_series_scale^(0.5), 'l')

## MSE is constant for fixed (b-a)^2, n
unif_series_shift <- vector("double", 100)
for (i in seq_along(unif_series_shift)) {
    unif_series_shift[i] <- unif_mle_mse(i, i+1, 10)
}
plot(1:100, unif_series_shift, 'l')



# Ex3
ex3_data <- c(3.23, -2.50, 1.88, -0.68, 4.43, 0.17,
              1.03, -0.07, -0.01, 0.76, 1.76, 3.18,
              0.33, -0.31, 0.30, -0.61, 1.52, 5.43,
              1.54, 2.28, 0.42, 2.33, -1.03, 4.00,
              0.39)

ex3_mu <- mean(ex3_data)
ex3_sigma <- sd(ex3_data)
ex3_n <- length(ex3_data)
ex3_mle_tau <- ex3_mu + 1.64 * ex3_sigma
(ex3_tau_se_delta <- ex3_sigma * sqrt((2 + 1.64^2) / (2*ex3_n)))

ex3_bootstrap_sample <- matrix(rnorm(n=ex3_n * 1000, mean=ex3_mu, sd=ex3_sigma), ncol=ex3_n)
ex3_bootstrap_values <- apply(ex3_bootstrap_sample, 1, function(x) quantile(x, 0.95, names=FALSE))
ex3_bootstrap_est <- mean(ex3_bootstrap_values)
(ex3_bootstrap_se <- sd(ex3_bootstrap_values))

## Convergence plot
plot(1:length(ex3_bootstrap_values), sqrt(accumulate(ex3_bootstrap_values, function(x, y) x + (y - ex3_bootstrap_est)^2 / 1:length(ex3_bootstrap_values))), 'l')


## Ex 7d)
boot1 <- apply(matrix(rbernoulli(200*1000, 160/200), 1000), 1, mean)
boot2 <- apply(matrix(rbernoulli(200*1000, 148/200), 1000), 1, mean)

boot_psi <- boot1 - boot2
sd(boot_psi)


## Ex 9
ex9 <- rnorm(100, 5, 1)

se_delta <- exp(mean(ex9)) / sqrt(length(ex9))

ex9_bootstrap <- apply(matrix(rnorm(100*1000, mean(ex9), 1), nrow=1000), 1, function(x) exp(mean(x)))

se_bootstrap <- sd(ex9_bootstrap)

bootstrap <- function(x, b, f) {
    tboot <- vector("double", b)
    for (i in seq(b)) {
        xsample <- sample(x, length(x), replace = TRUE)
        tboot[i] <- f(xsample)
    }
    tboot
}
ex9_nonparametric <- bootstrap(ex9, 1000, function(x) exp(mean(x)))


theta_hat <- apply(matrix(rnorm(100*1000, 5, 1), nrow=1000), 1, function(x) exp(mean(x)))

t <- seq(min(theta_hat), max(theta_hat), length.out = 1000)
ex9_normal <- dnorm(t, exp(mean(ex9)), se_delta)

# TODO: plot ex9_bootstrap, ex9_normal, ex9_nonparametric, theta_hat
library(ggplot2)
ggplot() +
    geom_histogram(data=data.frame(x=theta_hat), mapping = aes(x=x, y=..density..), fill="red") +
    geom_histogram(data=data.frame(x=ex9_bootstrap), mapping = aes(x=x, y=..density..), fill="blue", alpha=0.5)

## TODO: Ex9.10:
## Analytical: (n*x^(n-1)/theta)
## c.f parametric, non-parametric bootstrap
n <- 50
thetha <- 1
