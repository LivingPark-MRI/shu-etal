# Load required libraries
if (!require('e1071')) install.packages('e1071'); library('e1071')
if (!require('crayon')) install.packages('crayon'); library('crayon')
if (!require('caret')) install.packages('caret'); library('caret')

# Load dataset
dataset = read.csv("outputs/radiomics/paperFeatures.csv")
dataset = subset(dataset, select=-c(X))

# Normalize data
min_max_norm <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }

dataset[,-ncol(dataset)] <- as.data.frame(lapply(dataset[,-ncol(dataset)], min_max_norm))

write.csv(dataset, "outputs/radiomics/normalizedData.csv")

cat(green("Succesfully normalized top features, index available in: outputs/radiomics/normalizedData.csv"))