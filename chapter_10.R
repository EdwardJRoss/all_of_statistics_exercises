## My own example
## We coss atoin 100 times and get 43 heads
## Wald p-value is ~ 0.84
## This implies that about 84% of the time the values should be less extreme
trials <- 10000
trial_sample <- apply(matrix(rbinom(100*trials, 1, 0.5), nrow=trials), 1, mean)
mean(abs(trial_sample - 0.5) <= 0.07)
## We get simulated result is about 86-87%


## Ex 10.7

twain_3 <- c(0.225, 0.262, 0.217, 0.240, 0.230, 0.229, 0.235, 0.217)
snodgrass_3 <- c(0.209, 0.205, 0.196, 0.210, 0.202, 0.207, 0.224, 0.223, 0.220, 0.201)

## a) Wald test for equality of means

## Plug in mean
tlw_estimate <- mean(twain_3) - mean(snodgrass_3)

## Plug in estimator standard deviation
## Three letter word
tlw_se <- sqrt(var(twain_3) / length(twain_3) + var(snodgrass_3) / length(snodgrass_3))

tlw_wald <- abs(tlw_estimate) / tlw_se

tlw_pvalue <- 2*pnorm(- tlw_wald)

tlw_95_ci <- c(tlw_estimate - 2*tlw_se, tlw_estimate + 2*tlw_se)

## These indicate the null hypothesis that the means are equal is very unlikely


## b) Use permutation test to avoid the use of large sample methods.

permute_statistic <- function(x, y, f) {
    combined <- c(x, y)
    permuted <- sample(combined)
    f(permuted[1:length(x)], permuted[length(x) + 1:length(y)])
}

permutation_sample <- function(x, y, f, b) {
    observed <- f(x, y)
    answer <- 0
    for (i in 1:b) {
        if (permute_statistic(x, y, f) > observed) {
            answer <- answer + 1/b
        }
    }
    answer
}

tlw_permute_pvalue <- permutation_sample(twain_3, snodgrass_3, function(x, y) mean(x) - mean(y), 10000)

## Ex 10.10
ex10_df <- data.frame(week = c(-2, -1, 1, 2), chinese = c(55, 33, 70, 40), jewish = c(141, 145, 139, 161))

## Is there enough evidence to refute the proportions are from the same distribution?
## Use the mean difference

ex10_df$chinese_prop <- ex10_df$chinese / mean(ex10_df$chinese)
ex10_df$jewish_prop <- ex10_df$jewish / mean(ex10_df$jewish)

permutation_sample(ex10_df$chinese_prop, ex10_df$jewish_prop, function(x, y) mean(x) - mean(y), 1000)

## This is small enough that we could do an exact permutation analysis
## But the approximation shows p-value ~ 0.4-0.6. Not strong enough evidence to refute the null hypothesis.
## Note the degrees of freedom here (3) is very small.

## Ex 10.10
ex11_df <- data.frame(treatment = c('placebo', 'chloropromazine', 'dimenhydrinate', 'pentobarbital (100 mg)', 'pentobarbital (150 mg)'),
                      patients = c(80, 75, 85, 67, 86),
                      nausea = c(45, 26, 52, 35, 37))

# These are Bernoulli(p, n). MLE of p
ex11_df$mle <- ex11_df$nausea / ex11_df$patients

# Estimated standard error
ex11_df$se <- sqrt(ex11_df$mle * (1 - ex11_df$mle) / ex11_df$patients)

# Difference from placebo
ex11_df$odds_ratio <- ex11_df$mle / ex11_df$mle[ex11_df$treatment == 'placebo']

ex11_df$delta <- ex11_df$mle - ex11_df$mle[ex11_df$treatment == 'placebo']
ex11_df$delta_se <- sqrt(ex11_df$se^2 + ex11_df$se[ex11_df$treatment == 'placebo'] ^2)

ex11_test <- ex11_df[ex11_df$treatment != 'placebo',]

ex11_test$wald <- abs(ex11_test$delta) / ex11_test$se

ex11_test$signif95 <- ex11_test$wald > qnorm(1 - (1 - 0.95) / 2)
# So Chloropromazine and Pentobarbital (150mg) have a significantly different mean to Placebo (lower)

ex11_test$bonfe_signif95 <- ex11_test$wald / 4 > qnorm(1 - (1 - 0.95) / 2)
# None are significantly different by the Bonferonni method

ex11_test$pvalue <- 2*pnorm( - ex11_test$wald)

ex11_test <- ex11_test[order(ex11_test$pvalue),]

ex11_test$bh <- seq(1:nrow(ex11_test)) / nrow(ex11_test)

ex11_test$bh_signif95 <- ex11_test$pvalue < ex11_test$bh * 0.05
# So Chloropromazine and Pentobarbital (150mg)have a significantly different mean to Placebo (lower) by BH method


## Ex 10.12
lambda_0 <- 1
n <- 20
alpha <- 0.05

## Simulate Poisson(lambda_0) and perform Wald test

wald_value <- function(t, t_0) {
    mu <- mean(t)
    sigma <- sd(t) / sqrt(length(t))
    abs(mu - t_0)  / sigma
}

trials <- 10000
wald_trials <- apply(matrix(rpois(20 * trials, lambda_0), nrow = trials), 1, function (x) wald_value(x, lambda_0))
sum(wald_trials > qnorm(1 - alpha / 2)) / length(wald_trials)
# ~ 0.07
## The false positive rate is only about 7%, not 5%.
## This is because we're not really in the asymptotic region
