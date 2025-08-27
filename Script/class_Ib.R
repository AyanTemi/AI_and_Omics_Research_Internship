getwd()

setwd("C:/Users/Admin/Documents/AI_Omics_Internship_2025/Module_I")

dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plot")

data <- read.csv(file.choose())

# View data in spreadsheet format
View(data)

# Check structure of the dataset
str(data)

# Convert 'gender' column to factor
data$gender <- as.factor(data$gender)
str(data)

# Convert 'diagnosis' column to factor
data$diagnosis <- as.factor(data$diagnosis)
str(data)

# Convert 'height' column to factor
data$smoker <- as.factor(data$smoker)
str(data)

# Convert factor to numeric using ifelse statement (Yes = 1, No = 0)
data$smoker <- ifelse(data$smoker == "Yes", 1, 0)
class(data$smoker)

write.csv(data, file = "clean_data/patient_info_clean.csv")

# Save the entire R workspace
save.image(file = "full_workspace.RData")

