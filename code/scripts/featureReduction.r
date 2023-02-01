# Load library
if (!require('mRMRe')) install.packages('mRMRe'); library('mRMRe')
if (!require('crayon')) install.packages('crayon'); library('crayon')

# Top number of features (user input)
args <- commandArgs(trailingOnly = TRUE)
numberOfFeatures = args[1]
cohortType = args[2]

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
if (cohortType == "reference"){
    df <- read.csv("outputs/radiomics/Reference Cohort/features.csv", header = TRUE)
}
if (cohortType == "multiple"){
    df <- read.csv("outputs/radiomics/Multiple 3T Scanner Cohort/features.csv", header = TRUE)
} 
if (cohortType == "pdstate"){
    df <- read.csv("outputs/radiomics/PD-state Cohort/features.csv", header = TRUE)
}

# Set output index as numeric value
df[,ncol(df)] <- as.numeric(df[,ncol(df)])

df = subset(df, select = -c(X))

# Create dataset
f_data <- mRMR.data(data = data.frame(df))

# Perform classic mRMR
solution <- mRMR.classic(data=f_data, target_indices=ncol(df), feature_count = as.integer(numberOfFeatures))

# Get index of top 7 features
if (cohortType == "reference"){
    write.csv(solutions(solution), "outputs/radiomics/Reference Cohort/featureIndex.csv")
    cat(green("Succesfully computed top features, index available in: outputs/radiomics/Reference Cohort/featureIndex.csv"))
}
if (cohortType == "multiple"){
    write.csv(solutions(solution), "outputs/radiomics/Multiple 3T Scanner Cohort/featureIndex.csv")
    cat(green("Succesfully computed top features, index available in: outputs/radiomics/Multiple 3T Scanner Cohort/featureIndex.csv"))
}
if (cohortType == "pdstate"){
    write.csv(solutions(solution), "outputs/radiomics/PD-state Cohort/featureIndex.csv")
    cat(green("Succesfully computed top features, index available in: outputs/radiomics/PD-state-cohort/featureIndex.csv"))
}