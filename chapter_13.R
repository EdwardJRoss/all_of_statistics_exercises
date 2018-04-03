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
