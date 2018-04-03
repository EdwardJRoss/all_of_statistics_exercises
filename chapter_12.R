# Ex6. Compare the risk of the MLE and James Stein estimator. Try various values of n and various vectors theta. Summarize your results.

js_est <- function(x) {
    max((1 - (length(x) - 2)/(sum(x*x))), 0)*x
}

js_risk <- function(x, theta) {
    dtheta = js_est(x) - theta
    sum(dtheta^2)
}

mle_risk <- function(x, theta) {
    dtheta = x - theta
    sum(dtheta^2)
}

cf_risk <- function(theta, nsamples=1000) {
    n <- length(theta)
    x <- matrix(rnorm(nsamples*n, theta, 1), nrow=n)
    js <- apply(x, 2, js_est)
    js_risk <- mean(apply((js - theta)^2, 2, sum))
    mle_risk <- mean(apply((x - theta)^2, 2, sum))
    c(js_risk, mle_risk)
}


# This takes a while
s <- seq(-10, 10, 1)
t <- as.matrix(expand.grid(t1=s, t2=s, t3=s))
a <- apply(t, 1, cf_risk)
df <- data.frame(t1=t[,1], t2=t[,2], t3=t[,3], js=a[1,], mle=a[2,])

## It looks like it falls fast and is symmetric around the origin
library(ggplot2)
ggplot(df) + geom_point(mapping=aes(x=t1, y=t2, size=log(mle/js)))

ggplot(df) + geom_line(mapping=aes(x=sqrt(t1^2+t2^2+t3^2), y=(mle/js)))


# Let's look at n, ||theta||
# This will take a while
dg <- expand.grid(n=c(3, 4, 10, 50, 100), ntheta = seq(0, 1000, 1))
# Calculate fixed norm
grisk <- do.call(rbind, lapply(mapply(rep, dg$ntheta / dg$n, dg$n), cf_risk))
dg$js <- grisk[,1]
dg$mle <- grisk[,2]

# The improvement at 0 norm(theta) = 0 is roughly n times
# The improvement looks "roughly" logistic
# Approximately twice as good where norm(theta) = n
# The difference becomes marginal where norm(theta) = n^2
ggplot(dg) + geom_line(aes(x=ntheta, y=mle), colour="red") + geom_line(aes(x=ntheta, y=js), colour="blue") + facet_grid(n~., scale="free")

# So when there are a lot of parameters, and the parameter is close to zero, James-Stein estimator is *much* better.
