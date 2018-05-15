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


# ii) QDA

spam_qda <- MASS::qda(is_spam ~ ., spam)


spam_qda_pred <- predict(spam_qda, spam)$class


table(spam_qda_pred, cspam)


# 33% false negative, 5% false positive
#             cspam
#spam_qda_pred    0    1
#            0 2101   82
#            1  687 1731

# iii) Logistic Regression

spam_log <- glm(is_spam ~ ., 'binomial', spam)

spam_log_pred <- predict(spam_qda, spam)$class


table(spam_log_pred, cspam)

# Identital to qda
#             cspam
#spam_log_pred    0    1
#            0 2101   82
#            1  687 1731

# iv) Classification Tree

library(rpart)

spam_tree <- rpart(is_spam ~ ., spam)

spam_tree_pred <- predict(spam_tree, spam) > 0.5


# This will be rediculously overfit
table(spam_tree_pred, cspam)

#              cspam
#spam_tree_pred    0    1
#         FALSE 2653  296
#         TRUE   135 1517
