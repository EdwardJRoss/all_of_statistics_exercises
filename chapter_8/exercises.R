# Example 8.2
nerve_frame <- read.table('nerve.dat', sep='', na.strings='', fill=T)

# Flatten
nerve <- as.matrix(nerve_frame)
dim(nerve) <- nrow(nerve_frame) * ncol(nerve_frame)
nerve <- nerve[!is.na(nerve)]

skewness <- function(data, index=1:length(data)) {
    x <- data[index]
    mean(((x - mean(x))/sd(x))^3)
}

bootstrap <- function(x, b, f) {
    tboot <- vector("double", b)
    for (i in seq(b)) {
        xsample <- sample(x, length(x), replace = TRUE)
        tboot[i] <- f(xsample)
    }
    tboot
}

bootstrap_ci_normal <- function(x, b, f, alpha) {
    se <- sd(bootstrap(x, b, f))
    expected <- f(x)
    c(expected - se * qnorm(1 - alpha/2), expected + se * qnorm(1 - alpha/2))
}

bootstrap_ci_pivotal <- function(x, b, f, alpha) {
    boot_x <- bootstrap(x, b, f)
    expected <- f(x)
    c(2*expected - quantile(boot_x, 1-alpha/2), 2*expected - quantile(boot_x, alpha/2))
}

bootstrap_ci_percentile <- function(x, b, f, alpha) {
    boot_x <- bootstrap(x, b, f)
    c(quantile(boot_x, alpha/2), quantile(boot_x, 1-alpha/2))
}



bootstrap_ci_normal(nerve, 1000, skewness, .05)

bootstrap_ci_pivotal(nerve, 1000, skewness, .05)

bootstrap_ci_percentile(nerve, 1000, skewness, .05)

## It can be useful to show the convergence
nerve_bootstrap <- bootstrap(nerve, 5000, skewness)

hist(nerve_bootstrap)

nerve_running_sd <- sqrt(purrr::accumulate(nerve_bootstrap, function (x, y) x + y^2) / seq(1:length(b)))

## We could do this with an accumulator with a clever merge sort
nerve_running_quantile <- vector("double", length(nerve_bootstrap))
for (i in 1:length(nerve_bootstrap)) {
    nerve_running_quantile[i] <- quantile(nerve_bootstrap[1:i], 0.975)
}

plot(nerve_running_sd)

plot(nerve_running_quantile)

## Using the boot package
boot::boot.ci(boot::boot(data=nerve, R=1000, statistic=skewness),
              type=c("norm", "perc", "bca"))

## Exercise 8.6
lsat <- c(576, 635, 558,578,666,580,555,661,651,605,653,575,545,572,594)
gpa <- c(3.39,3.30,2.81,3.03,3.44,3.07,3.00,3.43,3.36,3.13,3.12,2.74,2.76,2.88,3.96)
boot <- 1000

lsat_gpa_cor <- cor(lsat, gpa)

lsat_gpa_boot <- vector("double", length(lsat))
for (i in 1:boot) {
    index <- sample(1:length(lsat), replace=TRUE)
    lsat_gpa_boot[i] <- cor(lsat[index], gpa[index])
}

hist(lsat_gpa_boot)

## 95% confidence intervals
lsat_gpa_ci_normal <- c(lsat_gpa_cor - 2*sd(lsat_gpa_boot), lsat_gpa_cor + 2*sd(lsat_gpa_boot))

lsat_gpa_ci_pivotal <- c(2 * lsat_gpa_cor - quantile(lsat_gpa_boot, 0.975), 2 * lsat_gpa_cor - quantile(lsat_gpa_boot,0.0275))

lsat_gpa_ci_percentile <- c( quantile(lsat_gpa_boot, 0.0275),  quantile(lsat_gpa_boot,0.975))

## Ex 7.2

n <- 100
ex2_normal_radius <- vector("numeric", n)
ex2_normal_centre <- vector("numeric", n)
for (i in 1:n) {
    y <- rnorm(50)
    x <- exp(y)
    ex2_normal <- bootstrap_ci_normal(x, 100, skewness, .05)
    ex2_normal_radius[i] <- (ex2_normal[2] - ex2_normal[1]) / 2
    ex2_normal_centre[i] <- (ex2_normal[2] + ex2_normal[1])/2
}

## TODO: commentary
(ex2_pivotal <- bootstrap_ci_pivotal(y, 1000, skewness, .05))
(ex2_percentile <- bootstrap_ci_percentile(y, 1000, skewness, .05))


## Exmple 8.3
x <- rt(25, 3)
T <- function(x) (quantile(x, 0.75) - quantile(x, 0.25)) / 1.34


(ex2_normal <- bootstrap_ci_normal(x, 1000, T, .05))
(ex2_pivotal <- bootstrap_ci_pivotal(x, 1000, T, .05))
(ex2_percentile <- bootstrap_ci_percentile(x, 1000, T, .05))

## Examining Ex 5

iterate_dbl <- function(f, x0, n) {
    x <- x0
    ans <- vector("numeric", n)
    for (i in 1:n) {
        x <- f(x)
        ans[i] <- x
    }
    ans
}

n <- 100
iter <- 1000
## For plotting convergence
## Note is numerically unstable
unconditional_bootstrap_series <- iterate_dbl(function(x) x + mean(sample(rnorm(n), n, replace=TRUE))^2, 0, iter) / 1:iter

## This should be 1
(unconditional_bootstrap_series[iter] * n) / (2 - 1/n)


## Ex 6
n <- 100
mu <- 5
x <- rnorm(n, mu)

sample_size <- 1000
sample <- vector("double",  sample_size)
for (i in 1:sample_size) {
    sample[i] <- exp(mean(rnorm(n, mu)))
}


## Using if X_i ~ Normal(mu, 1) then sum(X_i) ~ Normal(n * mu, n)
distribution <- function(x) sqrt(n / 2*pi) * exp(-n*(log(x) - mu)^2) / x

t <- seq(min(sample), max(sample), by=0.01)

library(ggplot2)
ggplot() +
    geom_histogram(data = data.frame(x=sample), mapping=aes(x=x, y=..density..), bins=90) +
    geom_line(data = data.frame(x=t, y=distribution(t)), mapping=aes(x=x, y=y))


## Exercise 7
ex7_sample <- runif(50)

theta <- max(ex7_sample)

## Actual distribution: F(t) = t^n/theta^n; f(t) = n * t^(n-1) / theta^(n)

ex7_dist <- function(x, n=50, theta=1) n*x^(n-1) / theta^n

ex7_bootstrap <- bootstrap(ex7_sample, 1000, max)

t = seq(0, 1, by=0.01)

nbins <- 90
bin_width <- (max(ex7_bootstrap) - min(ex7_bootstrap)) / nbins


library(ggplot2)
ggplot() +
    geom_histogram(data = data.frame(x=ex7_bootstrap), mapping=aes(x=x, y=..density..), bins=90) +
    geom_line(data = data.frame(x=t, y=ex7_dist(t)), mapping=aes(x=x, y=y))

## Exercise 8

## Estimating the square of the mean (sqm)


sqm <- function(x) mean(x)^2



## Analytic bootstrap variance
vboot <- function(x, f=sqm) {
    bootstraps <- do.call(expand.grid, rep(list(x), length(x)))
    tdist <- apply(bootstraps, 1, f)
    mean(tdist ^ 2) - mean(tdist) ^ 2
}

## Approximate bootstrap variance
vboot_approx <- function(x, b=10000, f=sqm) {
    var(bootstrap(x, b, f))
}

var_n <- function(x, k) {
    mean((x-mean(x))^k) / length(x) ^ (k-1)
}

soln <- function(x) {
    xbar <- mean(x)
    n <- length(x)
    var_n(x, 4) + (2 - 3 / n) * var_n(x, 2)^2 + 4 * xbar * var_n(x,3) + 4 * xbar^2 * var_n(x,2)
}


x <- c(0, 0, 1)
vboot_approx(x)
vboot(x)
soln(x)

## Yay it works!
