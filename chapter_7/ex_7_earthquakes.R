library(tibble)
library(dplyr)
library(ggplot2)

empirical_distribution <- function (x) {
    function (t) {
        sapply(t, function(t) sum(x <= t) / length(x))
    }}

confidence_interval <- function(n, alpha) sqrt(log(2/a) / (2*n))

plot_distribution <- function(x, sample, alpha) {
    eps <- confidence_interval(length(sample), alpha)

    df <- tibble(t=x) %>%
        mutate(fhat = empirical_distribution(sample)(t),
               fhat_low = fhat - eps,
               fhat_high = fhat + eps)

    ggplot(df) +
        geom_ribbon(mapping=aes(x=t, ymin=fhat - eps, ymax = fhat + eps),
                    fill = "red", alpha = 0.2) +
        geom_line(mapping=aes(x=t, y=fhat))
}

# Exercise 7.7
# Data source: http://www.stat.cmu.edu/~larry/all-of-statistics/=data/fijiquakes.dat
earthquakes <- read.table("fijiquakes.dat", sep = "" , header = T,
                          na.strings ="", stringsAsFactors= F)

magnitudes <- earthquakes$mag

hist(magnitudes)

plot_distribution(seq(min(magnitudes), max(magnitudes), length.out = 1000),
                  magnitudes,
                  0.05)


## Result of Ex 7.6
est_diff <- function(cdf, a, b) cdf(b) - cdf(a)
est_diff_se <- function(cdf, a, b, n) {
    sqrt((cdf(b) - cdf(a) + (cdf(a) + cdf(b))^2)/n)
}
## Confidence interval in the asymptotic (normal) approximation
est_diff_ci <- function(sample, a, b, alpha) {
    if (b < a) { stop("b must be greater than a") }
    cdf <- empirical_distribution(sample)
    est <- est_diff(cdf, a, b)
    eps <- qnorm(1 - alpha/2) * est_diff_se(cdf, a, b, length(sample))
    c(est - eps, est + eps)
}

est_diff_ci(magnitudes, 4.3, 4.9, 0.05)
## So this says that there is a 95% chance that between 45% and 61% of Earthquakes
## near Fiji are greater than magnitude 4.3 and at most magnitude 4.9

# This agrees with naive point estimate
sum(magnitudes > 4.3 & magnitudes <= 4.9) / length(magnitudes)


## Exercise 7.8
## Data from http://www.stat.cmu.edu/~larry/all-of-statistics/=data/faithful.dat
faithful <- read.table("faithful.dat", sep = "" , header = T,
                       skip=25,
                       na.strings ="", stringsAsFactors= F)
waiting <- faithful$waiting

## It is clearly bimodal
hist(waiting)

## Plug-in estimates
(wait_n <- length(waiting))
(wait_mean <- mean(waiting))
(wait_se <- sd(waiting) / sqrt(wait_n))
(wait_ci_90 <- c(wait_mean - qnorm(0.95) * wait_se, wait_mean + qnorm(0.95) * wait_se))
## So there's a 90% chance that the wait time is actually
## between 69.5 minutes and 72.3 minutes

(wait_median <- median(waiting))


## Ex 7.10
cloud <- read.table("cloud_seeding.dat", sep = "" , header = T,
                       skip=14,
                       na.strings ="", stringsAsFactors= F)

## Estimate the difference in mean precipitation.
## Thist is just plug-in
(cloud_diff <- mean(cloud$Seeded_Clouds) - mean(cloud$Unseeded_Clouds))
## By design of experiment the two random vectors are uncorrelated
(cloud_diff_se <- sqrt((var(cloud$Seeded_Clouds) + var(cloud$Unseeded_Clouds)) / nrow(cloud)))
(cloud_ci_95 <- c(cloud_diff - qnorm(0.975) * cloud_diff_se,
                  cloud_diff + qnorm(0.975) * cloud_diff_se))

## This misses zero. There is a 95% chance this captured the parameter
## In which case there is an advantage to seeding.
