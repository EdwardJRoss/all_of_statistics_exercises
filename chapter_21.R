# Smoothing with orthogonal functions

glass <- read.table('glass.dat')

# cosine basis: phi_j(x) = sqrt(2) cos(j * pi * x) on [0, 1]

# Estimation:
# Beta_j = 1/n sum_i phi_j(x)
# Risk: Rhat = sum(j = 1, J) sigma_j^2 / n + sum(j = J+1, sqrt(n)) (beta%2 - sigma_j^2 / n)+
# sigma_j^2 = 1/(n-1) * sum(i=1, n) (sigma_j(X_i) - beta_j)^2

# Regression:
# In the case of uniform points
# rhat(x) = sum_j=1^J beta_j phi_j(x)
# Risk: Rhat = J * sigma^2 / n + sum(j=J+1, n) * (b_j ^ 2 - simga^2 /n)+
# sigma%2 = n/k sum(i=n-k+1, n) beta_j^2, k = n/4
# Bands: c = RMS(phi_j(x)) * sigma * chi_J, a / sqrt(n)


# X = RI, Y = Al
# Exercise 7
# a) Do a nonparametric gregression to fit the modex Y = f(x) + eps
# The data are not on a regular grid, ignore this but do sort the data.
# Provide a function estimate, and estimate of the risk and a confidence bad.
# b) Use the wavelet method to estimate f

al <- glass$Al[order(glass$Al)]
ri <- glass$RI[order(glass$Al)]
almin <- min(al)
alrange <- diff(range(al))
# Normalise the x axis to the interval [0, 1]
al <- (al - almin) / alrange

cosine_basis <- function(j, x) {
    sqrt(2) * cos(j * pi * x)
}

# Estimate the coefficients
beta <- double(length(al))
for (i in seq_along(beta)) {
    beta[i] <- mean(ri*cosine_basis(i, al))
}

# Estimate the risk
n <- length(al)
k <- floor(n/4)
sigmahat2 <- n/k * sum(beta[seq(n - k + 1, n)])

risk <- function(j) {
    ans <- double(length(j))
    for (i in seq_along(j)) {
        varterm <- beta^2 - sigmahat2 / n
        var_est <- sum(varterm[seq(j[i]+1, n)])

        ans[i] <- j[i] * sigmahat2 / n + min(var_est, 0)
    }
    ans
}

all_risk <- risk(seq(1, floor(sqrt(length(al)))))

# Minimum at j = 6

jhat <- 6

reg_est <- function(x) {
    ans <- double(length(x))
    for (j in 1:jhat) {
        ans <- ans + beta[j] * cosine_basis(j, (x - almin) / alrange)
    }
    ans
}

# Our projected estimate
t <- seq(almin, almin + alrange, 0.01)
plot(glass$Al, glass$RI)
lines(t, reg_est(t))
rug(glass$Al)
rug(glass$RI, side=2)


# Error estimate?
a2 <- function(x) {
    ans <- double(length(x))
    for (j in 1:jhat) {
        ans <- ans + cosine_basis(j, (x - almin) / alrange)^2
    }
    ans
}

cbound <- function (x) {
    # 95% confidence interval
    qchisq(1 - 0.05, jhat) * sqrt(a2(x) * sigmahat2 / length(al))
}

library(dplyr)
library(ggplot2)

# I'm not sure I have implemented the bounds correctly, they seem very conservative...
glass %>% mutate(pred = reg_est(Al),
                 ymax=pred + cbound(Al),
                 ymin = pred - cbound(Al)) %>%
    ggplot(mapping = aes(x=Al, y = RI)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), fill="blue", alpha=0.3) +
    geom_point() +
    geom_rug() +
    geom_line(aes(y=pred))



# Haar wavelets:

# I'm getting numbers 1000 times too big!

# phi_j, k (x) = 2^(j/2) phi(2^j * x - k)
# phi(x) = {-1 0 <= x <= 1/2, 1 1/2 < x <= 1)

# Father wavelet
haarf <- function (x) { (x >= 0) * (x <= 1) }
# Mother wavelet
haarm <- function(x) { (x >= 0) * (x  <=1) * ifelse(x <= 0.5, -1, 1) }
# This is the general wavelet formula; with haarm replaced by mother wavelet
haar <- function(j, k, x) { 2^((j-1)/2) * haarm(2^(j-1) * x - (k-1)) }


haar_regression <- function(x, y) {
    n <- length(x)
    if (n != length(y)) {
        stop("x and y must be same length")
    }

    # Estimate the coefficients
    # Obviously haarf(al) == 1, but shows structure for different wavelet
    alpha <- mean(y * haarf(x))

    D <- list()
    for (i in 1:floor(log(n, 2))) {
        D[[i]] <- double(2^(i-1))
        for (j in seq_along(D[[i]])) {
            D[[i]][j] <- mean(y * haar(i, j, x))
        }
    }

    # It's not clear to me that this estimate is any good when the distribution is
    # far from normal (and if it's normal why are you using wavelets?)

    # Unfortunately this thresholds everything to zero!
    sigmahat <- sqrt(n) * median(abs(D[[i]])) / 0.6745

    #sigmahat <- 0

    # Use universal thresholding; this is similar to a BIC approach
    # Trying to create the smallest model instead of most accurate
    haar_beta <- list()
    for (i in 1:floor(log(n, 2))) {
        haar_beta[[i]] <- double(2^(i-1))
        for (j in seq_along(haar_beta[[i]])) {
            if (abs(D[[i]][j]) > sigmahat * sqrt(2*log(n) / n)) {
                haar_beta[[i]][j] <- D[[i]][j]
            }
        }
    }

    list(beta=haar_beta, alpha=alpha)
}

wavelet_reg <- function(haar_coeff, x) {
    ans <- haar_coeff$alpha * haarf(x)
    for (i in seq_along(haar_coeff$beta)) {
        for (j in seq_along(haar_coeff$beta[[i]])) {
            coeff <- haar_coeff$beta[[i]][j]
            ans <- ans + coeff * haar(i, j, x)
        }
    }
    ans
}

glass_haar <- haar_regression(al, ri)

# Doesn't look very good...
t <- seq(almin, almin + alrange, 0.01)
plot(glass$Al, glass$RI)
lines(t, wavelet_reg(glass_haar, (t - almin) / alrange ))
rug(glass$Al)
rug(glass$RI, side=2)



# Exercise 9 b

doppler <- function(x) { sqrt(x*(1-x)) * sin ( 2.1 * pi / (x + 0.05)) }

t <- seq(1/1024, 1, 1/1024)
y <- doppler(t) + 0.1 * rnorm(length(t))


doppler_haar <- haar_regression(t, y)


plot(t, y)
lines(t, wavelet_reg(doppler_haar, t))


# Exercise 11: Start
n <- 100
nsamples <- 1000
x <- matrix(rnorm(n * nsamples), nrow = n)

samples <- apply(x, 2, function(x) median(abs(x)))
mean(samples)
