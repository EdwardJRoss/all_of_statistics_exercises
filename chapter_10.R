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
    f(permuted[1:length(x)], permuted[length(x) + 1:length(combined)])
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
