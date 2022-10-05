# Load library
if (!require('mRMRe')) install.packages('mRMRe'); library('mRMRe')
if (!require('crayon')) install.packages('crayon'); library('crayon')

# Top number of features (user input)
numberOfFeatures <- commandArgs(trailingOnly = TRUE)

check.integer <- function(N){
    !length(grep("[^[:digit:]]", as.character(N)))
}

if (length(numberOfFeatures) > 1) { 
    cat(red("Error found! Too many arguments. Please supply 1 argument.\n"))
    stop() 
}
if (length(numberOfFeatures) == 0) { 
    cat(red("Error found! Missing argument. Please supply 1 argument.\n"))
    stop()
}
if (check.integer(numberOfFeatures) != TRUE) {
    cat(red("Error found! The argument must be an integer.\n"))
    stop()
}

# Read CSV
df <- read.csv("outputs/radiomics/features.csv", header = TRUE)

# Set output index as numeric value
df[,ncol(df)] <- as.numeric(df[,ncol(df)])

# Create dataset
f_data <- mRMR.data(data = data.frame(df))

# Perform classic mRMR
solution <- mRMR.classic(data=f_data, target_indices=ncol(df), feature_count = as.integer(numberOfFeatures))

# Get index of top 7 features
write.csv(solutions(solution), "outputs/radiomics/featureIndex.csv")

cat(green("Succesfully computed top features, index available in: outputs/radiomics/featureIndex.csv"))