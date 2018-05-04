# Non parametric curve estimation

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
        last_x = x[i]

        while (index <= length(breaks) && x[i] > breaks[index]) {
            if (!(last_break < breaks[index])) {
                stop(paste0('Break points not sorted at ', index))
            }
            last_break = breaks[index]

            index = index + 1
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
    qnorm(1/2 + (1/2 - alpha / 2)^(1/m)) / (2 * sqrt(n*h))
}

# Confidence intervals
## Phi(1/2 - (1/2-alpha/2)^1/m) / (2 sqrt(n*h))

# alpha / 2m

# (1 - x)^(1/m) = 1 - x/m

library(KernSmooth)

####

# 2 / (n-1) h - (n+1)/(n-1) * sum{j=1}^{m}(p_j ^2)


# KDE: 1/(h * n^2) * sum(i, j) K*(X_i - X_j / j) + 2/ (n*h)* K(0)
# K*(x) = K^(2) (x) - 2 * K(x)
# K^(2)(x) = int (K(z+y) * K(y)) dy


# f(x) = 1/n sum(1/h * K(x - X_i / h)
# K epanechinikov: 3/4 ( 1 - x^2/5) / sqrt(5) if |x| < sqrt(5)



# Exercise 3: Perform a nonparametric regression to fit RI = f(Al) + eps
# Use cross-validation to estimate the bandwidth
# Construct 95% CI
