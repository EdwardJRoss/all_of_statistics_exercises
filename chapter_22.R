spam_path = 'spam.data'
# See https://web.stanford.edu/~hastie/ElemStatLearn/datasets/spam.info.txt for description
if (!file.exists(spam_path)) {
    download.file('https://web.stanford.edu/~hastie/ElemStatLearn/datasets/spam.data',
                  spam_path)
}


spam <- read.csv(spam_path, sep=' ', header=FALSE)
colnames(spam)[ncol(spam)] <- 'is_spam'

mspam <- as.matrix(spam[,-ncol(spam)])
cspam <- spam[['is_spam']]


# Exercise 3a)


# MASS::qda takes a different approach
# dx = (x - mu) / sqrt(n)
# Do a QR decomposition dx = QR (Q unitary, R triangula)
# Then calculate: Scaling = (R )^{-1} by backpropogation (Solving triangular R * Scaling = I)
# Then S = (dx)^t (dx) =  R^T R
# So: x^T S^{-1} x' = (scaling x')^t * (scaling x)
# This is a really clever (and efficient) way to calculate an inner product with a different metric!
# (And of course it's more robust)
discriminant_analysis <- function(x, y) {
    n <- nrow(x)
    # Check length(y) == n
    p <- ncol(x)

    categories <- sort(unique(y))

    k <- length(categories)


    sigma <- array(dim = c(k, p, p))
    mu <- array(dim = c(k, p))
    pi <- double(k)
    for (i in 1:k) {
        xi <- x[y == categories[i],]
        pi[i] <- nrow(xi) / n
        mu[i, ] <- apply(xi, 2 ,mean)
        # Should probably use the MLE which differs by (n-1)/n
        sigma[i, , ] <- cov(xi)
    }

    list(categories = categories, pi = pi, sigma = sigma, mu = mu)
}

qda_predict <- function(d, x) {
    detsigma <- apply(d$sigma, 1, det)

    p <- ncol(x)
    # Could check p conforms with sigma/mu
    n <- nrow(x)


    k <- length(d$categories)

    delta0 <- log(d$pi) - 0.5 * log(detsigma)

    delta <- matrix(nrow=k, ncol=n)
    # Could convert to matrix operations (see predict.qda in MASS)
    # Loops can't be parallelized as well
    for (i in 1:k) {
        for (j in 1:n) {
            #delta[i, j] <- dx[,j] %*% sigmainv_dx[,j] + delta0[i]
            dx <- x[j,] - d$mu[i,]
            sigmainv_dx <- solve(d$sigma[i,1:p,1:p], dx)
            delta[i, j] <- delta0[i] - 0.5 * dx %*% sigmainv_dx
        }
    }

    pred <- apply(delta, 2, which.max)

    d$categories[pred]
}

lda_predict <- function(d, x) {
    # Average of correlation matrix
    # weighted by probability of outcome
    sigma <- apply(d$pi * d$sigma, c(2, 3), sum)


    k <- length(d$categories)
    n <- nrow(x)

    delta <- matrix(nrow=k, ncol=n)
    for (i in 1:k) {
        dx <- t(x) - 0.5 * d$mu[i,]
        # This is a constant, we could cache this
        sigmainv_mu <- solve(sigma, d$mu[i,])
        delta[i,] <- log(d$pi[i]) + t(dx) %*% sigmainv_mu
    }

    pred <- apply(delta, 2, which.max)

    d$categories[pred]
}


# A test set small enough to do QDA on by hand
test_data <- data.frame(x1 = c(0, 1, 0, 2, 1, 0, 1),
                        x2 = c(0, 0, 1, 0, 1, 2, -1),
                        y  = c(1, 1, 1, 0, 0, 0, 0))
x <- test_data[-3]
y <- test_data[[3]]
d <- discriminant_analysis(x, y)
md <- MASS::qda(x, as.factor(y))



dspam <- discriminant_analysis(mspam, cspam)


# i) LDA

spam_lda <- MASS::lda(is_spam ~ ., spam)

spam_lda_pred <- predict(spam_lda, spam)$class

# Confusion matrix
table(spam_lda_pred, cspam)

# Error around 11%

# 5% false negatives, 27% false positives

#            cspam
#spam_lda_pred    0    1
#            0 2663  387
#            1  125 1426

# By hand
spam_lda_pred_hand <- lda_predict(dspam, mspam)

table(spam_lda_pred, spam_lda_pred_hand)


# ii) QDA

spam_qda <- MASS::qda(is_spam ~ ., spam)


spam_qda_pred <- predict(spam_qda, spam)$class


table(spam_qda_pred, cspam)


# In some sense QDA actually does worse!
# 33% false negative, 5% false positive
#             cspam
#spam_qda_pred    0    1
#            0 2101   82
#            1  687 1731

# By hand
spam_qda_pred_hand <- qda_predict(dspam, mspam)
table(spam_qda_pred, spam_qda_pred_hand)


# iii) Logistic Regression

spam_log <- glm(is_spam ~ ., 'binomial', spam)

spam_log_pred <- predict(spam_qda, spam)$class

# Predictions are identical to QDA!
table(spam_log_pred, spam_qda_pred)

# Error is around 17%


# iv) Classification Tree

# rpart is only one possible library
# There's also e.g. "tree" which also uses the CART algorithms
# and C50 which uses Quinlan's C5.0 algorithm
library(rpart)

spam_tree <- rpart(is_spam ~ ., spam, method="class")


library(rpart.plot)
rpart.plot(spam_tree, type=4, extra=101)

nm <- colnames(predict(spam_tree))
spam_tree_pred <- apply(predict(spam_tree), 1, which.max)

# This will be rediculously overfit
table(nm[spam_tree_pred], cspam)

#              cspam
#spam_tree_pred    0    1
#         FALSE 2653  296
#         TRUE   135 1517


# We should really prune the tree using cross-validation

printcp(spam_tree)

prune(spam_tree, cp=0.03)

# Construct by hand?

## This finds the best split for a single predictor (column) and single categorical factor (outcome)
## using the Gini index
split_variable <- function(column, outcome) {
    idx <- order(column)
    lvls <- levels(outcome)

    unik <- rev(!duplicated(rev(column[idx])))

    run <- vapply(levels(outcome), function(x) x == outcome[idx], numeric(length(outcome)))

    runsum <- apply(run, 2, cumsum)
    totals <- apply(run, 2, sum)

    size <- sum(totals)

    gini <- function(row) {
        prop = row / sum(row)
        row_complement = totals - row
        prop_complement = row_complement / sum(row_complement)

        # This is the "prior" of a random varaible being in each split
        # which we use to weight the Gini index
        # This isn't described in all of statistics, but is used in the rpart alogorithm
        weight <- sum(row) / size

        # Gini for each split is (1 - sum(prob_k^2))
        # Total gini is the weighted sum of the two
        weight * (1 - sum(prop^2)) + (1 - weight) * (1 - sum(prop_complement^2))
    }

    gini_coeff <- apply(runsum[unik,], 1, gini)
    best_idx <- which.min(gini_coeff)
    best_gini <- gini_coeff[best_idx]

    # Calculate Gini, which.min
    # which to lookup unik
    best_split <- column[idx][which(unik)[best_idx]]

    one_length <- sum(runsum[unik,][best_idx,])
    two_length <- size - one_length
    length <- min(one_length, two_length)

    list(split=best_split, gini=best_gini, length=length)
}

best_split <- function(covariate, outcome) {
    splits <- apply(covariate, 2, function(x) split_variable(x, as.factor(outcome)))

    split_var <- which.min(lapply(splits, function(x) x$gini))
    ans <- splits[[split_var]]
    ans[["colnum"]] <- split_var
    ans
}

# This phase of the algorithm is a recursive procedure to
# split down the tree until some stopping condition is met
# (e.g. max depth, bucket size too small, "complexity" not decreasing enough)
split_root <- best_split(mspam, cspam)

left <- mspam[,split_root$colnum] <= split_root$split
right <- !left

splitl <- best_split(mspam[left,], cspam[left])
splitr <- best_split(mspam[right,], cspam[right])



# b) Apply 5-fold cross-validation

# This is built into modelr/rsample in the tidyverse
# as well as carat (and probably a dozen other libraries)
library(modelr)
library(dplyr)
library(purrr)

spam_folds <- crossv_kfold(spam, 5)

set.seed(126)
spam_folds %>% mutate(lda = map(train, function(df) MASS::lda(is_spam ~ ., df)))

pred_err <- function(pred, outcome) {
    matches <- pred != outcome
    mean(matches)
}

outcome <- function(model) {
    f <- model$call[['formula']]
    depvar <- all.vars(f)[[1]]
    model[[depvar]]
}

spamf <- as.formula(is_spam ~ .)

# When I use QDA I tend to get rank deficincy
# Maybe because some of the predictors are correlated?
spam_fits <- spam_folds %>%
    mutate(lda_mod = map(train, ~ MASS::lda(spamf, .)),
           log_mod = map(train, ~ glm(is_spam ~ ., 'binomial', .)),
           outcome = map(test, ~ as.data.frame(.)$is_spam),
           lda_pred = map2(lda_mod, test, ~ predict(.x, newdata=.y)$class),
           log_pred = map2(log_mod, test, ~ predict(.x, newdata=.y) > 0.5),
           lda_err = map2_dbl(lda_pred, outcome, pred_err),
           log_err = map2_dbl(log_pred, outcome, pred_err))

# Both around 10-11%
# For logistic regression it does very well *except* in fold 5
# it looks like there's something there that's bringing it down (correlated predictors?)
spam_fits %>% summarise(lda_sd = sd(lda_err),
                        log_sd = sd(log_err),
                        log_err = mean(log_err),
                        lda_err = mean(lda_err))

# By hand


# https://stackoverflow.com/questions/7659833/inverse-of-which
invwhich<-function(indices, outlength, useNames = TRUE)
{
    rv<-logical(outlength)
    if(length(indices) > 0)
    {
        rv[indices]<-TRUE
        if(useNames) names(rv)[indices]<-names(indices)
    }
    return(rv)
}

# We round n to the nearest folds
# Then return a list of indexes for each cross-fold
cv_mask <- function(n, folds=5) {
    size <- n %/% folds


    fold_seq <- seq(0, folds - 1)

    cv_idxs <- sample(seq(1, folds*size))

    lapply(fold_seq, function(x) invwhich(cv_idxs[seq(size * x, size * (x+1))], n))
}


set.seed(128)
cv_masks <- cv_mask(nrow(spam))

# This is pretty gnarly; probably clearer as a loop
# Outcome is fairly consistent
lda_error <- vapply(cv_masks, function(subset) mean(predict(MASS::lda(is_spam ~ ., spam[!subset,]), spam[subset,])$class != spam[subset, 'is_spam']), double(1))
log_error <- vapply(cv_masks, function(subset) mean( (predict(glm(is_spam ~ ., family='binomial', spam[!subset,]), spam[subset,]) > 0.5) != spam[subset, 'is_spam']), double(1))

# 11% for LDA and 8% for Logistic
mean(lda_error)
mean(log_error)


## Exercise 3c)
zero <- cspam == 0
one <- cspam == 1

# Calculate the p values for whether the means are unequal
spam_t <- apply(mspam, 2, function(x) t.test(x[zero], x[one])$p.value)

top10 <- order(spam_t) %in% seq(1,10)
lda10_error <- vapply(cv_masks, function(subset) mean(predict(MASS::lda(is_spam ~ ., spam[!subset,c(top10, TRUE)]), spam[subset,])$class != spam[subset, 'is_spam']), double(1))
log10_error <- vapply(cv_masks, function(subset) mean( (predict(glm(is_spam ~ ., family='binomial', spam[!subset,c(top10, TRUE)]), spam[subset,]) > 0.5) != spam[subset, 'is_spam']), double(1))
# It's actually made our predictors significantly worse

mean(lda10_error)
mean(log10_error)


# Let's try the strategy from Max Kuhn, Applied Predictive Modelling, Section 3.6
# Remove the most correlated predictors
# 1. Calcualte the correlation matrix of the predictors
# 2. Determine the two predictors associated with largest absolute pairwise correlation
# 3. For each determine the average correlation with other variables
# 4. Remove the variable with highest average correlation
# 5. Repeat until there are no absolute correlations above the threshold


threshold <- 0.75
removed <- c()

spamcor <- abs(cor(spam) - diag(ncol(spam)))

while (max(spamcor) > threshold) {
    cor_vars <- arrayInd(which.max(spamcor), dim(spamcor))

    vara <- cor_vars[,1]
    varb <- cor_vars[,2]

    cora <- mean(spamcor[vara,])
    corb <- mean(spamcor[varb,])

    if (corb > cora) {
        v <- varb
    } else {
        v <- vara
    }

    removed <- c(removed, v)
    spamcor <- spamcor[-v, -v]
}

kept_vars <- !invwhich(removed, ncol(mspam))

ldacor_error <- vapply(cv_masks, function(subset) mean(predict(MASS::lda(is_spam ~ ., spam[!subset,c(kept_vars, TRUE)]), spam[subset,])$class != spam[subset, 'is_spam']), double(1))
logcor_error <- vapply(cv_masks, function(subset) mean( (predict(glm(is_spam ~ ., family='binomial', spam[!subset,c(kept_vars, TRUE)]), spam[subset,]) > 0.5) != spam[subset, 'is_spam']), double(1))

# About the same or slightly worse
mean(ldacor_error)
mean(logcor_error)


# 5. Classify Spam Data using SVMs
# We could alternatively use implementation at http://svmlight.joachims.org/
# There is an R interface to this implementation in the klaR package

# install.packages('e1071')
svmlinear_spam <- e1071::svm(is_spam ~ ., spam, type='C-classification', kernel='linear')
svmlinear_pred <- predict(svmlinear_spam, spam)

table(cspam, svmlinear_pred)
#     svmlinear_pred
#cspam    0    1
#    0 2663  125
#    1  185 1628

# I'm aware there are some built in validation tools in SVM, but for fairness use k-fold CV to evaluate
# Could make this a bit shorter with a function
library(e1071)

svm_fits <- spam_folds %>%
    mutate(linear_mod = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='linear')),
           outcome = map(test, ~ as.data.frame(.)$is_spam),
           linear_pred = map2(linear_mod, test, ~ predict(.x, newdata=as.data.frame(.y))),
           linear_err = map2_dbl(linear_pred, outcome, pred_err))

# A linear model actually does slightly better than logistic regression
svm_fits %>% summarise_at(vars(ends_with("err")), mean)


# Excercise 22.6

# LDA is a linear classifier so VC dimension is d+1
d <- ncol(iris) - 1
n <- nrow(iris)

# Error is > 1
# Too loose a bound to be informative!
alpha <- 0.05
vc_error95 <- sqrt(32 / n * log(8 * (n^(d+1) + 1) / alpha))

iris_lda <- MASS::lda(Species ~ ., iris)
training_loss <- mean(predict(iris_lda)$class != iris$Species)

# Excercise 22.8

# Could probably do some clever mapping here?
svmpoly_fits <- spam_folds %>%
    mutate(mod1 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=1, gamma=1, coef0=1)),
           mod2 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=2, gamma=1, coef0=1)),
           mod3 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=3, gamma=1, coef0=1)),
           mod4 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=4, gamma=1, coef0=1)),
           mod5 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=5, gamma=1, coef0=1)),
           mod6 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=6, gamma=1, coef0=1)),
           mod7 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=7, gamma=1, coef0=1)),
           mod8 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=8, gamma=1, coef0=1)),
           mod9 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=9, gamma=1, coef0=1)),
           mod10 = map(train, ~ svm(spamf, ., type = 'C-classification', kernel='polynomial', degree=10, gamma=1, coef0=1)),
           outcome = map(test, ~ as.data.frame(.)$is_spam))

svmpoly_preds <- svmpoly_fits %>%
    mutate_at(vars(starts_with('mod')), .funs=funs(map2(., test, ~ predict(.x, newdata=as.data.frame(.y))))) %>%
    select(starts_with('mod'), 'outcome')

svmpoly_cverr <- mutate_at(svmpoly_preds, vars(starts_with('mod')), .funs=funs(map2_dbl(., outcome, pred_err))) %>%
    select(-outcome)


summarise_all(svmpoly_cverr, mean)

summarise_all(svmpoly_cverr, sd)
