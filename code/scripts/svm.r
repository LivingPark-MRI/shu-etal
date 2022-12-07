# Load required libraries
if (!require('e1071')) install.packages('e1071'); library('e1071')
if (!require('randomForest')) install.packages('randomForest'); library('randomForest')
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

# Fit to SVM
classifier = svm(formula = output ~ .,
                 data = trainSet,
                 type = 'C-classification',
                 kernel = 'linear')

y_pred <- predict(classifier, testSet[,-ncol(testSet)])
roc_object <- roc( as.numeric(testSet[,ncol(testSet)]), as.numeric(y_pred))
auc(roc_object)