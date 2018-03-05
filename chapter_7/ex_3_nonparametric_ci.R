library(tibble)
library(dplyr)
library(ggplot2)

empirical_distribution <- function (x) {
    function (t) {
        sapply(t, function(t) sum(x < t) / length(x))
    }}

confidence_interval <- function(n, alpha) sqrt(log(2/alpha) / (2*n))

bound_interval <- function(sample, cdf, alpha) {
    fhat <- empirical_distribution(sample)
    eps <- confidence_interval(length(sample), alpha)
    all(abs(fhat(sample) - cdf(sample)) < eps)
}

plot_distribution <- function(x, sample, cdf, alpha) {
    eps <- confidence_interval(length(sample), alpha)

    df <- tibble(t=x) %>%
        mutate(f = cdf(t),
               fhat = empirical_distribution(sample)(t),
               fhat_low = fhat - eps,
               fhat_high = fhat + eps)

    ggplot(df) +
        geom_ribbon(mapping=aes(x=t, ymin=fhat - eps, ymax = fhat + eps),
                    fill = "red", alpha = 0.2) +
        geom_line(mapping=aes(x=t, y=f)) +
        geom_line(mapping=aes(x=t, y=fhat))
}

# How often does the alpha confidence interval contain n-estimator of CDF?
run_trials <- function(generator, cdf, n, alpha, trials) {
    cdf <- pnorm
    generator <- rnorm

    results <- vector("logical", trials)
    for (j in seq(trials)) {
        results[j] <- bound_interval(generator(n), cdf, alpha)
    }
    sum(results) / length(results)
}

run_trials(rnorm, pnorm, 100, 0.05, 1000)
# 0.97

run_trials(rcauchy, pcauchy, 100, 0.05, 1000)
# 0.96

plot_distribution(seq(-2.5, 2.5, by=0.01), rcauchy(100), pcauchy, 0.05)
