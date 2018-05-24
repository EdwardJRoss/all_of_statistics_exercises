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


# iv) Classification Tree

library(rpart)

spam_tree <- rpart(is_spam ~ ., spam)


library(rpart.plot)
rpart.plot(spam_tree, type=4)

spam_tree_pred <- predict(spam_tree, spam) > 0.5


# This will be rediculously overfit
table(spam_tree_pred, cspam)

#              cspam
#spam_tree_pred    0    1
#         FALSE 2653  296
#         TRUE   135 1517

# Construct by hand?
# We need to check for each variable, and for each possible split value
# the Gini impurity (which is


# b) Apply 5-fold cross-validation

n <- nrow(spam)
set.seed(428)
cv_idxs <- shuffle(1:n)

function cv(model, sample) {

}
