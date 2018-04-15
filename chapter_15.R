# Exercise 5
# Age and financial status
# Age: 1 = <35, 2= 35-54, 3 = 55+
# Financial status: 1 = worse, 2= same, 3 = better than a year ago
df5 <- read.csv('montanadat.dat', sep='\t', skip=19, na='*')

# Base r solution
t5 <- table(df5$AGE, df5$FIN)

# At a glance it seems older people were worse off
# Let's check whether they are actually associated

# Chisquared test in base r
chisq.test(t5)

loglike_indep <- function(t) {
    m1 <- margin.table(t, 1)
    m2 <- margin.table(t, 2)
    total <- margin.table(t)

    ans <- 0
    for (i in seq(nrow(t))) {
        for (j in seq(ncol(t))) {
            ans <- ans + 2 * t[i, j] * log(t[i,j] * total / (m1[i] * m2[j]))
        }
    }
    ans
}

log5 <- loglike_indep(t5)
pval5 <- 1 - pchisq(log5, nrow(t5) + ncol(t5) - 2)

# They are correlated with p-value 2e-4 (1 in 5000)

# Using a library
DescTools::GTest(t5)

# dplyr solution
library(dplyr)
library(tidyr)
indep_test <- df5 %>%
    rename_all(tolower) %>%
    filter(!is.na(age), !is.na(fin)) %>%
    count(age, fin) %>%
    group_by(age) %>%
    mutate(margin_age = sum(n)) %>%
    group_by(fin) %>%
    mutate(margin_fin = sum(n)) %>%
    ungroup() %>%
    mutate(total = sum(n),
           expected = margin_age * margin_fin / total,
           loglik = 2 * n * log(n * total / (margin_fin * margin_age)),
           chi2 = (n - expected)^2 / expected) %>%
    summarise(loglik = sum(loglik), chi2 = sum(chi2),
              nage = n_distinct(age), nfin = n_distinct(fin)) %>%
    gather(test, value, loglik, chi2) %>%
    mutate(pval = 1 - pchisq(value, nage + nfin - 2))

# All indicators agree: They are associated!

# Exercise 6
# Estimate correlation beween temperature and latitude

# Of course if there is a non-linear relationship correlation isn't best measure
df6 <- read.csv('USTemperatures.dat', sep='\t', skip=15)

# Roll your own method
# Non-parametric plug in estimate (With bias correction)
# (Bias correction should be about 2% difference for 50 observations)
lat <- df6$Lat
sdlat <- sqrt(sum((lat - mean(lat))^2) / (length(lat) - 1))
temp <- df6$JanTemp
sdtemp <- sqrt(sum((temp - mean(temp))^2) / (length(temp) - 1))
corr <- sum((lat - mean(lat)) * (temp - mean(temp))) / (sdlat * sdtemp * (length(lat) - 1))

# Built in
cor(lat, temp)

# Pearson correlation confidence interval
pearson_corr_ci <- function(corr, n, alpha = 0.05) {
    fcorr <- 0.5 * (log(1 + corr) - log(1 - corr))

    sehat <- qnorm( 1 - alpha / 2) / sqrt(n - 3)

    fci <- c(fcorr - sehat, fcorr + sehat)

    (exp(2*fci) - 1) / (exp(2 * fci) + 1)


}

# Strong negative correlation
pearson_corr_ci(corr, length(lat))

# Base r
cor.test(temp, lat)

# Bootstrap confidence interval
mycor <- function(df, x, y, indices) {
    cor(df[indices, x], df[indices, y])
}

boot6 <- boot::boot(data=df6, statistic=mycor, x="Lat", y="JanTemp", R=1000)

# All in the range -0.66 to -0.99
boot::boot.ci(boot6, conf = 0.95)

# Exercise 7
# Test whether calcium intake and drop in blood pressure are associated
df7 <- read.csv('Calcium.dat', sep='\t', skip=16)


# Use the two Sample Kolmogorov-Smirnov test

# From first principles
ca_decrease <- df7$Decrease[df7$Treatment == 'Calcium']
pl_decrease <- df7$Decrease[df7$Treatment == 'Placebo']

# To compare the Empircal Distribution Functions we just need
# to know how many points of one are to the left of the other
nca <- length(ca_decrease)
npl <- length(pl_decrease)

edf <- function(data) {
    function(x) {
        unlist(lapply(x, function(x) sum(data <= x)/ length(data)))
    }
}

ca_edf <- edf(ca_decrease)
pl_edf <- edf(pl_decrease)

jump_points <- sort(unique(c(ca_decrease, pl_decrease)))

# D statistic
# This could be calculated a little more efficiently by skipping the edf
d6 <- max(abs(ca_edf(jump_points) - pl_edf(jump_points)))

# Numerically unstable
H <- function(t, eps = 1e-6) {
    ans <- double(length(t))
    for (i in seq_along(t)) {
        j <- 1
        ans[i] <- 1
        finished <- FALSE
        while (!finished) {
            term <- 2 * (-1)^(j) * exp(- 2 * j^2 * t[i]^2)
            j <- j + 1
            if (abs(term) < eps) {
                finished <- TRUE
            }
            ans[i] <- ans[i] + term
        }
    }
    ans
}

# Solve by bisection
bisect_solve <- function(f, a, b, eps = 1e-6) {
    mid <- (a + b) / 2
    fl <- f(a)
    fr <- f(b)
    if (fl * fr > 0)
        stop("Bad signs")
    fmid <- f(mid)

    if (abs(fr - fl) > eps) {
        if (fl * fmid < 0)
            bisect_solve(f, a, mid, eps)
        else
            bisect_solve(f, mid, b, eps)
    } else {
        mid
    }
}
# Test: bisect_solve(function(x) x*x - 2, 0, 2)

# I'm wrong here!

h95 <- bisect_solve(function(x) H(x) - 0.95, 1, 10)

sqrt(nca * npl / (nca + npl)) * d6
h95

# We don't have enough evidence to reject null hypothesis with 95% confidence
1 - H(sqrt(nca * npl / (nca + npl)) * d6)
# p-value: 0.345

# Built-in R
ks.test(ca_decrease, pl_decrease)
