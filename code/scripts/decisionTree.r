# Load required libraries
if (!require('e1071')) install.packages('e1071'); library('e1071')
if (!require('rpart.plot')) install.packages('rpart.plot'); library('rpart.plot')
if (!require('crayon')) install.packages('crayon'); library('crayon')
if (!require('caret')) install.packages('caret'); library('caret')
if (!require('pROC')) install.packages('pROC'); library('pROC')

# Load dataset
dataset = read.csv("outputs/radiomics/normalizedData.csv")

# Split train/test sets
sample <- sample(nrow(dataset), size=floor(.70*nrow(dataset)), replace=F)
trainSet <- dataset[sample, ]
testSet  <- dataset[-sample, ]

# Remove useless columns
trainSet = subset(trainSet, select=-c(X))
testSet = subset(testSet, select=-c(X))

# Classify
classifier <- rpart(output ~ ., data = trainSet)

y_pred <- predict(classifier, testSet[,-ncol(testSet)])
roc_object <- roc( as.numeric(testSet[,ncol(testSet)]), as.numeric(y_pred))
print(auc(roc_object))