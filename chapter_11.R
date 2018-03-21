## Simulation (11.4)
param_range <- c(0,1)
prior <- function(x) 1
lhood <- function(x, p) dbinom(x, 1, p)


posterior <- function(x, data) {
    ans <- vector("double", length(x))
    for (i in seq_along(x)) {
        ans[i] <- prior(x[i]) * prod(lhood(data, x[i]))
    }
    ans
}


data <- c(1, 0, 0, 1, 0, 0, 1)

est_mean <- integrate(function (x) x*posterior(x, data), 0, 1)$value / integrate(function (x) posterior(x, data), 0, 1)$value

exact_mean <- (sum(data) + 1) / (length(data) + 2)

phi <- function(x) log(x / (1-x))

est_phi_mean <- integrate(function(x) phi(x) * posterior(x, data), 0, 1)$value /  integrate(function (x) posterior(x, data), 0, 1)$value


## This disagrees with Wasserman's answer in Example 11.3 because he dropped an exp(phi) from the numerator of the derivative.
exact_phi_dist <- function(x) dbeta(exp(x)/(1+exp(x)), sum(data) + 2, length(data) - sum(data) + 2)

exact_phi_mean <- integrate(function(x) x * exact_phi_dist(x), -100, 100)$value / integrate(exact_phi_dist, -100, 100)$value


## Ex 11.2
ex11_2 <- rnorm(100, 5, 1)

## Take f(mu) = 1 and find the posterior density.
ex11_2_posterior <- function(mu) {
    n <- length(ex11_2)
    xbar <- mean(ex11_2)
    sqrt(n/(2*pi)) * exp(-n/2*(mu - xbar)^2)
}

t <- seq(0, 10, by=0.05)
plot(t, ex11_2_posterior(t), 'l')

# Similar
ex11_2_simulation <- rnorm(1000, mean(ex11_2), 1/sqrt(length(ex11_2)))
hist(ex11_2_simulation)

# d) Find posterior density for exp(mu)

ex11_2_exp_density <- function(t) {
    ans <- vector("double", length(t))
    for (i in seq_along(t)) {
        ans[i] <- 1/ t[i] * prod(1/ex11_2 * exp(-(exp(ex11_2) - t[i])^2/ 2))
    }
    ans
}
