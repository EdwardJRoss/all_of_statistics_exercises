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
        ans[i] <- 1/ t[i] * prod(exp(-(ex11_2 - log(t[i]))^2 / 2))
    }
    ans
}

## Simulation
hist(exp(ex11_2_simulation))

# e-f) Test against simulation

ex2_mu <- 5
ex2_n <- 100
ex2_sample_size <- 10000

ex2_samples <- matrix(rnorm(ex2_sample_size * ex2_n, ex2_mu, 1), ncol = ex2_n)

# posterior interval
ex2_in_pi <- function(x, mu=ex2_mu, alpha=0.05) {
    centre <- mean(x)
    size <- qnorm(1 - alpha/2) / sqrt(length(x))
    # c(centre - size, centre+size)
    abs(centre - mu) < size
}


## 95.16%: This seems excellent
sum(apply(ex2_samples, 1, ex2_in_pi)) / ex2_sample_size

ex2_in_exp_ci <- function(x, mu=ex2_mu, alpha=0.05) {
    centre <- exp(mean(x))
    size <- qnorm(1 - alpha/2) * exp(mean(x)) / sqrt(length(x))
    # c(centre - size, centre+size)
    abs(centre - exp(mu)) < size
}

# 94.85%, acceptable as is approximate confidence interval
sum(apply(ex2_samples, 1, ex2_in_exp_ci)) / ex2_sample_size


# Ex 4b)
n1 <- 50
s1 <- 30
n2 <- 50
s2 <- 40

p1 <- s1/n1
p2 <- s2/n2

boot_size <- 1000
bootstrap1 <- rbinom(boot_size, n1, p1) / n1
bootstrap2 <- rbinom(boot_size, n2, p2) / n2

bootstrap_se <- sd(bootstrap2 - bootstrap1)

bootstrap_90ci <- c((p2 - p1) - qnorm(1 - 0.1/2) * bootstrap_se, (p2 - p1) + qnorm(1 - 0.1/2) * bootstrap_se)


# c) Use prior f(p1, p2) = 1. Use simulation to find posterior mean and 90% Posterior interval.

sim_size <- 1000
# Prior for each of p1 and p2 is Beta(1, 1) and they are indep
# Posterior is separately Beta(s + 1, n - s + 1)
sim_post1 <- rbeta(sim_size, s1 + 1, n1 - s1 + 1)
sim_post2 <- rbeta(sim_size, s2 + 1, n2 - s2 + 1)

tau_post <- sim_post2 - sim_post1
hist(tau_post)

# A little lower, but similar to the bootstrap
post_90ci <- quantile(tau_post, c(0.05, 0.95))


# d) Check: se = sqrt(1/(n1*p1*(1-p1)) + 1/(n2*p2*(1-p2)))
sqrt(1/(n1*p1*(1-p1)) + 1/(n2*p2*(1-p2)))
sd(log((bootstrap1)/(1-bootstrap1) / ((bootstrap2) / (1-bootstrap2))))

# e)

phi_post <- log((sim_post1) / (1-sim_post1) / ((sim_post2) / (1-sim_post2)))

# Bayesian mean
mean_phi  <- mean(phi_post)

# Very close to Frequentist approach
post_90ci_phi <- quantile(phi_post, c(0.05, 0.95))

# Ex 5
n <- 10
s <- 2


t <- seq(0, 1, 0.01)

alpha <- c(0.5, 1, 10, 100)
beta <- c(0.5, 1, 10, 100)

ex5_df <- data.frame(
    t = rep(t, length(alpha)),
    alpha = rep(alpha, each = length(t)),
    beta = rep(beta, each = length(t))
    )

ex5_df$posterior <- with(ex5_df, dbeta(t, alpha + s, beta + n - s))
ex5_df$alpha_beta <- with(ex5_df, paste(as.character(alpha), as.character(beta), sep=","))

# The stronger the prior the more evidence that is needed to overcome it
library(ggplot2)
ggplot(ex5_df, mapping=aes(x=t, y=posterior, group = alpha_beta, col=alpha_beta)) + geom_line()
