# Non parametric curve estimation

# See http://vita.had.co.nz/papers/density-estimation.pdf
# for an overview of tools in R


# Exercise 2: Estimate density of Refractive Index using histogram and kernel density estimator.
# Use cross-validation to choose the amount of smoothing.
# Experiment with different binwidths and bandwidths.
# Construch 95% confidence bands.

glass <- read.table('glass.dat')

# breaks is ordered vector of *right* endpoints
# Always adds a default endpoint of +infinity
# x is data points in sorted order
# Returns an integer vector one greater than length of breaks
# with the number of x between those two points
histogram_vector_explicit <- function(x, breaks) {
    ans <- integer(length(breaks) + 1)
    index <- 1
    last_break <- -Inf
    last_x <- -Inf
    for (i in seq_along(x)) {

        if (!(last_x <= x[i])) {
            stop(paste0('Data points not sorted at ', i))
        }
        last_x <- x[i]

        while (index <= length(breaks) && x[i] > breaks[index]) {
            if (!(last_break < breaks[index])) {
                stop(paste0('Break points not sorted at ', index))
            }
            last_break <- breaks[index]

            index <- index + 1
        }
        ans[index] <- ans[index] + 1
    }
    ans
}

# A simpler implementation taking advantage of R
# Differs from above implementation by not adding infinities
# and breaks can be an integer for number of slices to cut into
histogram_vector <- function(x, breaks) {
    table(cut(x, breaks))
}

# Explicit cross-validation risk as a check
histogram_risk_cv <- function(x, m) {
    cuts <- cut(x, m)
    h <- diff(range(x)) / m
    n <- length(x)
    fhat <- table(cuts) / (n * h)

    v <- sum(fhat^2 * h)

    fhati <- double(length(x))
    for (i in seq_along(x)) {
        # fhat calculated excluding Xi, evaluated at Xi
        fhati[i] <- sum(cuts[-i] == cuts[i]) / ((n-1) * h)
    }
    b <- sum(fhati)

    v - (2/n) * b
}


# Cross validation risk for length 1 m
histogram_risk1 <- function(x, m) {
    n <- length(x)
    h <- diff(range(x)) / m
    prob <- histogram_vector(x, m) / n

    2 / (h*(n-1)) - sum(prob^2) * (n+1)/(h*(n-1))
}

# Estimate risk of histogram estimator of x with m pieces
# using cross validation
# m can be a vector
histogram_risk <- function(x, m) {
    n <- length(x)
    h <- diff(range(x)) / m
    counts <- lapply(m, function(m) histogram_vector(x, m))

    # Prob = count / n
    prob2 <- vapply(counts, function(x) sum((x/n)^2), double(1))

    (2  - prob2 * (n+1)) / (h * (n-1))
}

RI <- glass$RI
l <- seq(2, length(RI))
r <- histogram_risk(RI, l)

# It is fairly flat around the minimum
plot(l, r, 'l')

best_bin <- l[which.min(r)]
oversmoothed <- round(best_bin / 5)
undersmoothed <- round(best_bin * 5)

hist(RI, best_bin)

hist(RI, undersmoothed)

hist(RI, oversmoothed)

hist_ci <- function(x, m, alpha) {
    h <- diff(range(x)) / m
    n <- length(x)
    qnorm(0.5 + 0.5 * (1 - alpha)^(1/m)) / (2 * sqrt(n*h))
}

library(ggplot2)
library(dplyr)

# Plotting a histogram from scratch using ggplot
rihist <- hist(RI, breaks=best_bin, plot=FALSE)
ridata <- with(rihist,
               data.frame(left=breaks[-length(breaks)],
                          right=breaks[-1],
                          counts=counts,
                          density=density)) %>%
    mutate(upper95 = density + hist_ci(RI, best_bin, 0.05),
           lower95 = vapply(ridata$lower95, function(x) max(x, 0), double(1)))

ggplot(ridata) + geom_rect(aes(xmin=left, xmax=right, ymin=0, ymax=density))

ggplot(ridata) +
    geom_rect(aes(xmin=left, xmax=right, ymin=lower95, ymax=upper95)) +
    geom_segment(aes(x=left, xend=right, y=density, yend=density), colour="black") +
    geom_segment(aes(x=left, xend=left, y=lag(density), yend=density), colour="black")

# It would be interesting to know how sharp this calculated 95% CI is

# An alternative approach would be to use a (semiparametric?) bootstrap
# We treat each bin like an independent Binomial distribution
# and resample the histogram

# Both these approaches assume the bins are constant (hence the "semiparametric?")
# This is ok here since the minimum is so flat.


# Kernel approach

# Base R
plot(density(RI, bw="ucv", kernel="epanechnikov"))

plot(density(RI, bw="ucv", kernel="gaussian"))

# First principles
# Epanechinkov kernel
epan <- function(x) {
    ifelse(abs(x) < sqrt(5), 3/(4*sqrt(5)) * (1 - x^2/5), 0)
}

plot(t, epan(t), 'l')

# Returning a function probably isn't the best decision
# A better choice would be to return an S3 object
# that contained the relevant infomation:
# * The points (x)
# * The weights (1/(n*h))
# * The scale h
# * The kernel
# We could also potentially optimise this case by not adding
# points outside the bandwidth, making prediction more efficient

kern_dens <- function(x, h, kernel=epan) {
    n <- length(x)
    function (t) {
        ans <- double(length(t))
        for (i in seq_along(t)) {
            ans[i] <- sum(kernel((t[i] - x) / h)) / (n*h)
        }
        ans
    }
}

t <- seq(-5, 5, by=0.01)
plot(t, kern_dens(RI, 0.3, kernel=epan)(t), 'l')

plot(t, kern_dens(RI, 0.3, kernel=dnorm)(t), 'l')

# The direct but inefficient way
# Jhat = integral(fhat^2) - 2/n * sum_i (fhat_{-i} (x_i))
kern_cv_risk_direct <- function(x, hs, kernel=dnorm) {
    n <- length(x)
    ans <- double(length(hs))
    for (i in seq_along(hs)) {
        h <- hs[i]

        fit2 <- function(t) (kern_dens(x, h, kernel)(t))^2

        term1 <- integrate(fit2, min(x), max(x), subdivisions = 1000)$value

        term2 <- 0
        for (j in seq_along(x)) {
            term2 <- term2 + kern_dens(x[-j], h, kernel)(x[j])
        }
        ans[i] <- term1 - (2 / n) * term2
    }
    ans
}

t <- seq(0.01, 2, 0.01)
cvrisk <- kern_cv_risk_direct(RI, t)
plot(t, cvrisk, 'l')

hmin <- t[which.min(cvrisk)]


# Cross validation estimate of the risk
# Use Gaussian kernel for simplicity
# This doesn't work, need to investigate why
broken_gauss_kern_cv_risk <- function(x, h) {
    n <- length(x)
    # K(x) convolved with K(-x)
    k2 <- function (x) dnorm(x, sd=2)
    kstar <- function(x) k2(x) - 2*dnorm(x)

    cterm <- 2 * dnorm(0) / (n*h)

    ans <- double(length(h))
    for (i in seq_along(h)) {
        for (xi in x) {
            ans[i] = ans[i] + sum(kstar((xi - x) / h[i])) / (h[i] * n^2)
        }
    }
    ans + cterm
}

# TODO: Extend to general kernels using FFT for efficiency
# (Note kernels must integrate to 1, and have 0 mean, and positive standard deviation)
# Re(fft(y*y, inverse=TRUE))


# Oversmoothed
t <- seq(-5, 5, by=0.01)
plot(t, kern_dens(RI, hmin * 6, kernel=dnorm)(t), 'l')

# Undersmoothed
t <- seq(-5, 5, by=0.01)
plot(t, kern_dens(RI, hmin / 6, kernel=dnorm)(t), 'l')

# Just right
t <- seq(-5, 5, by=0.01)
plot(t, kern_dens(RI, hmin, kernel=dnorm)(t), 'l')


# There are other ways to do this too:
# The key observation that the ideal h is proportional to n^(-4/5) means if we have a good
# estimate for a sample we can extend that to the population
# KernSmooth takes an approach combining binning and "Direct plug in estimator"
# Kernel Smoothing - Wand and Jones (1995)
# There's also the "Sheather-Jones" method used in density (bandwidth="SJ") (also in the above book)
# Sheather, S. J. and Jones, M. C. (1991) A reliable data-based bandwidth selection method for kernel density estimation. Journal of the Royal Statistical Society series B, 53, 683–690.
# https://www.jstor.org/stable/2345597
# See also: https://www.umiacs.umd.edu/labs/cvl/pirl/vikas/Software/optimal_bw/optimal_bw_code.htm

# For Confidence Intervals we could use the bootstap
# https://stats.stackexchange.com/questions/207129/compute-confidence-interval-for-univariate-kernel-density-estimator
# https://cran.r-project.org/web/packages/HSAUR/vignettes/Ch_density_estimation.pdf
# The formula in the book is heuristic anyway

# However we could also do it via estimating each point as a Binomial trial

# It's easier to use the vectorised representation in density than mine as a closure

B <- 2500
bw <- 0.24

kernd <- density(RI, n=length(RI), bw=bw, kernel="epanech")
fromx <- min(kernd$x)
tox <- max(kernd$x)

estimates <- matrix(NA, nrow=length(RI), ncol=B)
for (b in 1:B) {
    xstar <- sample(RI, replace=TRUE)
    dstar <- density(xstar, n=length(RI), bw=bw, kernel="epanech", from=fromx, to=tox)
    estimates[,b] <- dstar$y
}

ConfidenceBands <- apply(estimates, 1, quantile, probs = c(0.025, 0.975))

kerndf <- data.frame(x = kernd$x,
                     y = kernd$y,
                     ymin = ConfidenceBands[1,],
                     ymax = ConfidenceBands[2,])



library(ggplot2)

ggplot(kerndf, aes(x=x, y=y)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), fill="blue", alpha=0.3) +
    geom_line() +
    geom_rug(aes(x=RI, y=as.double(NA)), data=glass)

# Of course this estimate is very optimistic where there is no data

# The famous KernSmooth package
library(KernSmooth)


plot(bkde(RI, bandwidth=bw), type='l')

# The 'oversmoothed bandwidth selector' seems to do a bad job
plot(bkde(RI), type='l')


# The direct plugin seems to do an ok job
plot(bkde(RI, bandwidth=dpik(RI)), type='l')

# It also predicts a binwidth that coalesces into a single peak
hist(RI, round(diff(range(RI)) / dpih(RI)))

# Exercise 3: Perform a nonparametric regression to fit RI = f(Al) + eps
# Use cross-validation to estimate the bandwidth
# Construct 95% CI


# Base R
# It is probably only valid near the centre where points are close together
plot(ksmooth(glass$Al, glass$RI, kernel="normal", bandwidth=0.3), type='l')
points(glass$Al, glass$RI)

# We'll estimate the risk using cross validation
# This is really only valid if the points are pretty even
# Remove outliers

al <- glass$Al[glass$Al <= 3 & glass$Al >= 0.7]
ri <- glass$RI[glass$Al <= 3 & glass$Al >= 0.7]


# I'll do this the "pedestrian" way rather than using the fast shortcut formula
# This doesn't work: The output error looks flat
broken_kreg_cv_risk <- function(x, y, h, kernel="normal") {
    if (length(x) != length(y)) { stop("x and y must be same length") }
    risk <- double(length(h))

    for (i in seq_along(h)) {
        for (j in seq_along(x)) {
            smooth <- ksmooth(x[-j], y[-j], bandwidth=h[i], kernel=kernel, range.x = range(x),
                              n.points = length(x))
            risk[i] <- risk[i] + (y[j] - smooth$y[j])^2
        }
    }
    risk
}

t <- seq(.3, 6, .3)
krisk <- kreg_cv_risk(al, ri, t)

plot(t, krisk)

bw = 0.3

# We'll use the bootstrap to estimate the error
ri_al_smooth <- ksmooth(al, ri, kernel="normal", range.x = range(al), n.points=length(al), bandwidth=bw)

B = 1000
estimates <- matrix(NA, nrow=length(al), ncol=B)
for (b in 1:B) {
    indices <- sample(seq(1, length(al)), replace=TRUE)
    xstar <- al[indices]
    ystar <- ri[indices]
    dstar <- ksmooth(xstar, ystar, kernel="normal", range.x = range(al), n.points=length(al), bandwidth=bw)
    estimates[,b] <- dstar$y
}

# The na.rm is a little dodgy here
# When we take a bootstrap sample with the x-axis too bunched up we may not be able to predict
# value in region with no x, so returns NA
# Does this introduce a bias? Maybe.
ConfidenceBands <- apply(estimates, 1, quantile, probs = c(0.025, 0.975), na.rm=TRUE)

kreg_df <- data.frame(x = ri_al_smooth$x,
                     y = ri_al_smooth$y,
                     ymin = ConfidenceBands[1,],
                     ymax = ConfidenceBands[2,])


ggplot(kreg_df, aes(x=x, y=y)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), fill="blue", alpha=0.3) +
    geom_line() +
    geom_rug(aes(x=Al, y=RI), data=glass)
