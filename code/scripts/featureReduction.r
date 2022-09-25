# Load library
library(mRMRe)

# Read CSV
df <- read.csv("outputs/radiomics/features.csv", header = TRUE)

# Set output index as numeric value
df[,ncol(df)] <- as.numeric(df[,ncol(df)])

# Create dataset
f_data <- mRMR.data(data = data.frame(df))

# Perform classic mRMR
solution <- mRMR.classic(data=f_data, target_indices=57, feature_count = 7)

# Get index of top 7 features
write.csv(solutions(solution), "outputs/radiomics/featureIndex.csv")

print("Succesfully computed top features, index available in: outputs/radiomics/featureIndex.csv")