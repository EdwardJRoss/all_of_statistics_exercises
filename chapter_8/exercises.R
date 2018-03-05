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

nerve_running_sd <- vector("double", length(nerve_bootstrap))
nerve_running_quantile <- vector("double", length(nerve_bootstrap))
for (i in 1:length(nerve_bootstrap)) {
    nerve_running_sd[i] <- sd(nerve_bootstrap[1:i])
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

lsat_gpa <- matrix(c(lsat, gpa), ncol=2)

lsat_gpa_corr <- corr(lsat_gpa)

lsat_gpa_boot <- vector("double", nrow(lsat_gpa_boot))
for (i in 1:boot) {
    index <- sample(1:nrow(lsat_gpa_boot))
    lsat_gpa_boot[i] <- corr(lsat_gpa[index, ])
}

bootstrap_ci_normal(lsat_gpa, 1000, corr, .05)

bootstrap_ci_pivotal(nerve, 1000, skewness, .05)

bootstrap_ci_percentile(nerve, 1000, skewness, .05)
