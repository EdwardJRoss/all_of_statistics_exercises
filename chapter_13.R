# From http://lib.stat.cmu.edu/DASL/Datafiles/carmpgdat.html
#    VOL: Cubic feet of cab space
#    HP: Engine horspower
#    MPG: Average miles per gallon
#    SP: Top speed (mph)
#    WT: Vehicle weight (100 lb)
carmpg <- read.csv('carmpg.csv')


linear_regression <- function(x, y) {
    solve(t(x) %*% x, t(x) %*% y)
}

## Exercide 6a)
## Fit linear regressoin MPG from HP.

mod_6a <- lm(MPG ~ HP, carmpg)

# The higher the HP the lower the MPG; makes sense
est_6a <- linear_regression(cbind(1, carmpg$HP), carmpg$MPG)


# carmpg$pred_hp <- predict(mod_6a1, carmpg)
carmpg$pred_6a <- as.vector(cbind(1, carmpg$HP) %*% est_6a)

# carmpg$resid_6a <- residuals(mod_6a)
carmpg$resid_6a <- carmpg$MPG - carmpg$pred_6a

# Estimate of variance
var_6a <- sum(carmpg$resid_6a^2) / (length(carmpg$MPG) - 2)

se_6a0 <- sqrt(var_6a * sum(carmpg$MPG^2)) / (length(carmpg$MPG) * sd(carmpg$MPG))
se_6a1 <- sqrt(var_6a / length(carmpg$MPG)) / sd(carmpg$MPG)

se_6a <- c(se_6a0, se_6a1)

coeff_6a_bound <- cbind(est_6a - se_6a, est_6a + se_6a)

wald_6a = est_6a / se_6a


estvar_6a <- function(x) {
    n <- length(carmpg$HP)
    factor <- var_6a / (n * var(carmpg$HP))
    ans <- vector("double", length(x))
    for (i in seq_along(x)) {
        ans[i] <- factor * (x[i]^2 - 2 * x[i] * mean(carmpg$HP) + sum(carmpg$HP^2)/n)
    }
    ans
}

carmpg$estvar_6a <- estvar_6a(carmpg$HP)

# Estimate of uncertainty in prediction value
predvar_6a <- function(x) {
    ans <- vector("double", length(x))
    for (i in seq_along(x)) {
        ans[i] <- var_6a * (1 + sum((carmpg$HP - x[i])^2) / (n * var(carmpg$HP)))
    }
    ans
}

library(ggplot2)
ggplot(carmpg) +
    # 95% confidence interval for points
    geom_ribbon(aes(x=HP, ymin=pred_6a - sqrt(estvar_6a), ymax = pred_6a + 2 * sqrt(estvar_6a)), alpha = 0.5) +
    geom_line(aes(x=HP, y=pred_6a)) +
    geom_point(aes(HP, MPG))

# Residuals
r2_6a <- sum(carmpg$resid_6a^2)

# The residuals clearly increase with increasing MPG; this makes the bounding analysis unlikely
ggplot(carmpg) + geom_point(aes(x=MPG, y=resid_6a))


# Ex 6b

mod_6b <- lm(MPG ~ log(HP), carmpg)

summary(mod_6b)

carmpg$pred_6b <- predict(mod_6b, carmpg)

carmpg$resid_6b <- residuals(mod_6b)


# Residuals are much smaller
r2_6b <- sum(carmpg$resid_6b^2)


# The fit looks a little better for large HP
ggplot(carmpg) +
    geom_line(aes(x=HP, y=pred_6a), colour="blue") +
    geom_line(aes(x=HP, y=pred_6b), colour="red") +
    geom_point(aes(HP, MPG))


# The residuals are closer, but still not even
ggplot(carmpg) +
    geom_point(aes(x=MPG, y=resid_6a), colour="blue") +
    geom_point(aes(x=MPG, y=resid_6b), colour="red")

# Canonical R with Confidence Interval

carmpg_6b <- cbind(carmpg, predict(mod_6b, interval="confidence", level=0.95))
ggplot(carmpg_6b, aes(x=HP, y=MPG, ymin=lwr, ymax=upr)) +
    geom_ribbon(alpha=0.1) +
    geom_line(aes(y=fit))  +
    geom_point()
# Exercise 7
carmpg_7 <- read.csv("carmpg.csv")

# First principles

linear_fit <- function(x, y) {
    solve(t(x) %*% x, t(x) %*% y)
}

# Doesn't make sense to use MAKEMODEL as it's sparse; coupld split out make
x7a <- as.matrix(cbind(INTERCEPT=1, carmpg_7[,c(-1,-4)]))
y7a <- carmpg_7[,4]

coeffs_7a <- linear_fit(x7a, y7a)
pred_7a <- x7a %*% coeffs_7a
resid_7a <- y7a - pred_7a

var_7a <- sum(resid_7a^2) / (nrow(x7a) - ncol(x7a))

# Sum of squares fit: 1000
r2_7a <- sum(resid_7a^2)

stderr_7a <- sqrt(var_7a * diag(solve(t(x7a) %*% x7a)))
# This is t value in r output
wald_7a <- coeffs_7a / stderr_7a

# This differs slightly from what R returns
# It seems to be using a slightly different test
pvalue_7a <- 2* pnorm(-abs(wald_7a))

# Spread looks ok
plot(x7a[,2], resid_7a)

plot(x7a[,3], resid_7a)

plot(x7a[,4], resid_7a)

# Using R
mod_7a <- lm(MPG ~ VOL + HP + SP + WT, carmpg_7)
summary(mod_7a)

plot(mod_7a)

# 7b)
# First principles
mallow_selector <- function(s) {
    x <- x7a[,s]
    fit <- linear_fit(x, y7a)
    pred <- x %*% fit
    resid <- y7a - pred
    rms <- sum(resid^2)
    if (is.null(nrow(x))) {
        c <- 1
        r <- length(x)
    } else {
        c <- ncol(x)
        r <- nrow(x)
    }
    var <- rms / (r - c)
    rms + 2 * length(s) * var
}

mallow_selector(seq(1,5))

# Backward Search
# Dropping 2 gives smallest result: 1060
mallow_selector(-1)
mallow_selector(-2)
mallow_selector(-3)
mallow_selector(-4)
mallow_selector(-5)

# These all increase
mallow_selector(c(-2, -1))
mallow_selector(c(-2, -3))
mallow_selector(c(-2, -4))
mallow_selector(c(-2, -5))

# So: MPG ~ 1 + HP + SP + WT


# Forward Search

# Smallest for 1: 8300
mallow_selector(1)
mallow_selector(2)
mallow_selector(3)
mallow_selector(4)
mallow_selector(5)

# Smallest for 1,5: 1500
mallow_selector(c(1, 2))
mallow_selector(c(1, 3))
mallow_selector(c(1, 4))
mallow_selector(c(1, 5))

# Smallest for 1, 5, 4: 1490.8
mallow_selector(c(1, 5, 2))
mallow_selector(c(1, 5, 3))
mallow_selector(c(1, 5, 4))


# 1, 5, 4, 3: 1141
mallow_selector(c(1, 5, 4, 3))
mallow_selector(c(1, 5, 4, 2))

# Increases
mallow_selector(c(1, 5, 4, 3, 2))

# So: MPG ~ 1 + WT + SP + HP
# Same result

# R automatic way
step(mod_7a, director="backward")
step(mod_7a, director="forward")


# Ex 7c): Zheng-Loh
order_7c <- order(-abs(wald_7a))

zl <- double(length(order_7c))
for (i in seq_along(order_7c)) {
    x <- x7a[,order_7c[seq(i)]]
    fit <- linear_fit(x, y7a)
    pred <- x %*% fit
    resid <- y7a - pred
    if (is.null(nrow(x))) {
        n <- length(x)
        c <- 1
    } else {
        n <- nrow(x)
        c <- ncol(x)
    }
    rms <- sum(resid^2)
    var <- rms / (n - c)
    zl[i] <- rms + i * var * log(n)
}
# Excludes 2: Same result


# d)

subsets <- function(x) {
    if (!is.null(x) && length(x) > 0) {
        ans <- list(x)

        for (i in seq_along(x)) {
            ans <- c(ans, subsets(x[-i]))
        }
        ans
    }
}

# This is a global maximum for Mallow
mallow_7d <- mapply(mallow_selector, subsets(seq(1, 5)))

colnames(x7a)[subsets(seq(1,5))[which.min(mallow_7d)][[1]]]


# We could always use basis expansion to get lots more variables!
# (logs interaction terms, etc)

linear_aic <- function(x, y, b) {
    n <- nrow(x)
    k <- ncol(x)
    rms <- sum((y - x %*% b)^2)

    n * log(2*pi) + n * log(rms / n) + 1 + 2*k
    #eps <- norm(y - x %*% b)

    #n/2 * log(n-k) - n * log(eps) - (n+k)/2
}

linear_bic <- function(x, y, b) {
    n <- nrow(x)
    k <- ncol(x)
    eps <- norm(y - x %*% b)

    n/2 * log(n-k) + k/2*log(n) - n * log(eps) - (n-k)/2
}

b <- linear_fit(x7a, y7a)
## TODO: BIC

bic_selection <- function(s) {
    x <- as.matrix(x7a[,s])
    b <- linear_fit(x, y7a)
    linear_bic(x, y7a, b)
}

bic_7d <- mapply(bic_selection, subsets(seq(1, 5)))
colnames(x7a)[subsets(seq(1,5))[which.min(bic_7d)][[1]]]

## BIC chooses a much more parsimonious model: HP
## AIC chooses: 1 + HP + SP + WT

## Exercise 11
coris <- read.csv('coris.dat', skip=3, header=FALSE,
                  col.names=c('rownum', 'sbp', 'tobacco','ldl',
                              'adiposity','famhist','typea','obesity',
                              'alcohol','age','chd'))[, -1]


## Reweighted Least Squares
## AKA Newton's Method
# Should probably monitor convergence
logistic_fit <- function(x, y, tol=1e-6, max_iter=1000) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    # Start with beta = 0. Since this function is pretty nice shoudl converge quickly.
    beta <- matrix(rep(0, ncol(x)))


    for (i in seq(max_iter)) {
        xbeta <- x %*% beta
        exp_b <- exp(xbeta)
        p <- exp_b / (1 + exp_b)

        z <- xbeta + (y - p) / (p * (1-p))

        w <- diag(as.vector(p * (1-p)))

        last_beta <- beta
        beta <- solve(t(x) %*% w %*% x, t(x) %*% w %*% z)

        sqerr <- sum((last_beta - beta)^2)
        if (sqerr < tol) {
            break
        }
    }

    if (i == max_iter) {
        stop(paste0("Too many iterations: ", sqerr))
    }
    beta
}

## Base R:
## mod <- glm(chd ~ ., family=binomial, coris)
## AIC(mod)
## MASS::stepAIC(mod)

x <- cbind(1, coris[-ncol(coris)])
y <- coris[ncol(coris)]
logistic_fit(x, y)

logistic_aic <- function(x, y, b) {
    x <- as.matrix(x)
    xb <- x %*% b
    loglik <- t(y) %*% xb - sum(log(1+exp(xb)))
    2 * (ncol(x) - loglik)
}

# Backwards search by AIC
aic_logistic_search <- function(x, y) {
    b <- logistic_fit(x, y)
    this_aic <- logistic_aic(x, y, b)
    aics <- double(ncol(x))
    for (i in seq(to=ncol(x))) {
        aics[i] <- logistic_aic(x[,-i], y, logistic_fit(x[,-i], y))
    }
    min_aic <- min(aics)
    min_index <- which.min(aics)

    if (this_aic <= min_aic) {
        b
    } else {
        # Optimisation: pass through b and this_aic
        aic_logistic_search(x[,-min_index], y)
    }
}

# Much slower than base R/MASS
best <- aic_logistic_search(x, y)

xb_best <- as.matrix(x[,rownames(best)]) %*% best

p_best <- exp(xb_best) / (1 + exp(xb_best))


best_error <- sum(abs(y - (p_best > 0.5))) / nrow(y)

# The best choice by backwards stepwise AIC (488) is:
# chd ~ tobacco + ldl + famhist + typea + age
# This is the most effective predictor; not implying causality
# The training error is 26% is quite high!
# Especially given p = 1 gives 35% error rate.


## Check on Exercise 13.9

num_sample_9 <- 10000
n_9 <- 1000
theta_9 <- 0

sample_9 <- matrix(rnorm(num_sample_9 * n_9, theta_9),
                   nrow = num_sample_9)

delta_9 <- apply(sample_9, 1, function(x) 1 - n/2 * mean(x)^2)

mean(delta_9 > 0)
