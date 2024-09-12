
#for EE vs EPL 

## Kw test with BH correction

library(stats)

T_Chemicals_data_EE_EPL2 <- read.csv("", check.names = FALSE)

# for numeric variables

T_Chemicals_data_EE_EPL2_numberic <- T_Chemicals_data_EE_EPL2[,-c(2,319,320)]

# Create an empty list for storing Kruskal-Wallis and BH correction results
results_list <- list()

# Create a vector to store variable names
variable_names <- colnames(T_Chemicals_data_EE_EPL2_numberic)[3:328]

for (i in 3:328) {
  # Execute Kruskal-Wallis test
  kw_result <- kruskal.test(T_Chemicals_data_EE_EPL2_numberic[, i], as.factor(T_Chemicals_data_EE_EPL2_numberic[, 2]))
  
  # Store Kruskal-Wallis test result in the list along with variable name
  result_entry <- list(Variable = variable_names[i - 2], Kruskal_Test = kw_result)
  results_list[[i - 2]] <- result_entry
}

# Extract p-values from the Kruskal-Wallis test results
p_values <- sapply(results_list, function(result) result$Kruskal_Test$p.value)

# Apply BH correction to the p-values
adjusted_p_values <- p.adjust(p_values, method = "BH")


# Add the adjusted p-values to the results list

for (i in 1:length(results_list)) {
  results_list[[i]]$adjusted_p_value <- adjusted_p_values[i]
}

# Convert the list to a data frame for easier writing to file
results_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(Variable = x$Variable, P_Value = x$Kruskal_Test$p.value, Adjusted_P_Value = as.numeric(x$adjusted_p_value))
}))

# Write the data frame to a CSV file
write.csv(results_df, "kw_with_BH_diagnosis_EE_EPLadditions_numberic.csv", row.names = FALSE)


# Filter the data frame to keep only features with adjusted p-value < 0.05
significant_results_dfEE_EPL_numberic  <- subset(results_df, Adjusted_P_Value < 0.05)
# Extract the names of significant features
significant_featuresEE_EPL_numberic <- significant_results_dfEE_EPL_numberic$Variable


columns_to_keep <- c("Stage", significant_featuresEE_EPL_numberic )


filtered_T_Chemicals_dataEE_EPL_numberic <- T_Chemicals_data_EE_EPL2_numberic[, significant_featuresEE_EPL_numberic, drop = FALSE]

# Select these columns from the transposed T_Chemicals_data
# Make sure "Stage" is a column in T_Chemicals_data
if("Stage" %in% colnames(T_Chemicals_data_EE_EPL2_numberic)) {
  filtered_T_Chemicals_dataEE_EPL_numberic <- T_Chemicals_data_EE_EPL2_numberic[, columns_to_keep, drop = FALSE]
} else {
  stop("Column 'Stage' not found in T_Chemicals_data_EE_EPL2_numberic")
}


# for category variables: The following results indicates no significant difference


# Explicitly perform Chi-square tests for each categorical variable

for (col_index in 319:320) {
  # Convert the column to a factor
  T_Chemicals_data_EE_EPL2[, col_index] <- factor(T_Chemicals_data_EE_EPL2[, col_index])
  
  # Ensure that the comparison column is also a factor if it is categorical
  comparison_column <- factor(T_Chemicals_data_EE_EPL2[, 3])  # Replace '2' with the correct column index if needed
  
  # Perform the Chi-square test
  chi_square_result <- tryCatch({
    chisq.test(table(T_Chemicals_data_EE_EPL2[, col_index], comparison_column))
  }, warning = function(w) {
    w
  }, error = function(e) {
    e
  }, finally = {
    cat("Processed column", col_index, "\n")
  })
  
  # Check if result is a chi-square test object and print p-value
  if (inherits(chi_square_result, "htest")) {
    cat("Variable:", colnames(T_Chemicals_data_EE_EPL2)[col_index], "- Chi-square p-value:", chi_square_result$p.value, "\n")
    # Add the Chi-square result to the results list
    results_list[[length(results_list) + 1]] <- list(
      Variable = colnames(T_Chemicals_data_EE_EPL2)[col_index],
      Chi_Square_Test = chi_square_result
    )
  } else {
    cat("Error or warning occurred in Chi-square test for column", col_index, "\n")
  }
}

# Perform Fisher's Exact Test for "location"

fisher_test_result <- tryCatch({
  fisher.test(table(T_Chemicals_data_EE_EPL2[, 319], T_Chemicals_data_EE_EPL2[, 2]))
}, warning = function(w) {
  cat("Warning in Fisher's test for column 319:", w$message, "\n")
}, error = function(e) {
  cat("Error in Fisher's test for column 319:", e$message, "\n")
})

# Print the p-value from Fisher's test
if (exists("fisher_test_result")) {
  cat("Variable:", colnames(T_Chemicals_data_EE_EPL2)[319], "- Fisher's Exact Test p-value:", fisher_test_result$p.value, "\n")
}

##------------------------------------------------
## EE vs EPL


library(dplyr) 
library(data.table) 
library(caTools)
library(pROC) 
library(ggplot2) 
library(ggpubr) 
library(ggprism) 
library(MASS)
library(caret) 

EE_EPL_addition = filtered_T_Chemicals_dataEE_EPL_numberic
EE_EPL_addition <-  EE_EPL_addition%>% mutate(Stage = ifelse(Stage =="EE",1,0))
EE_EPL_addition$Stage<-factor(EE_EPL_addition$Stage,levels = c(0,1),labels = c("EPL","EE"))

model_EE_EPL_addition <- glm(Stage ~ .,EE_EPL_addition , family = binomial)
model.both_EE_EPL_addition <- stepAIC(model_EE_EPL_addition, direction = 'both') # 'forward', 'backward'
# 查看模型结果
summary(model.both_EE_EPL_addition)$coefficients

##---------------------------------------------------------------------------------------
# adding demographic, and habitatary behabior ect. together with chemicals features

#for HC vs EPL 

## Kw test with BH correction

library(stats)

T_Chemicals_data_HC_EPL2 <- read.csv("", check.names = FALSE)

# for numeric variables

T_Chemicals_data_HC_EPL2_numberic <- T_Chemicals_data_HC_EPL2[,-c(2,319,320)]

# Create an empty list for storing Kruskal-Wallis and BH correction results
results_list <- list()

# Create a vector to store variable names
variable_names <- colnames(T_Chemicals_data_HC_EPL2_numberic)[3:328]

for (i in 3:328) {
  # Execute Kruskal-Wallis test
  kw_result <- kruskal.test(T_Chemicals_data_HC_EPL2_numberic[, i], as.factor(T_Chemicals_data_HC_EPL2_numberic[, 2]))
  
  # Store Kruskal-Wallis test result in the list along with variable name
  result_entry <- list(Variable = variable_names[i - 2], Kruskal_Test = kw_result)
  results_list[[i - 2]] <- result_entry
}

# Extract p-values from the Kruskal-Wallis test results
p_values <- sapply(results_list, function(result) result$Kruskal_Test$p.value)

# Apply BH correction to the p-values
adjusted_p_values <- p.adjust(p_values, method = "BH")


# Add the adjusted p-values to the results list

for (i in 1:length(results_list)) {
  results_list[[i]]$adjusted_p_value <- adjusted_p_values[i]
}

# Convert the list to a data frame for easier writing to file
results_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(Variable = x$Variable, P_Value = x$Kruskal_Test$p.value, Adjusted_P_Value = as.numeric(x$adjusted_p_value))
}))

# Write the data frame to a CSV file
write.csv(results_df, "kw_with_BH_diagnosis_HC_EPLadditions_numberic.csv", row.names = FALSE)


# Filter the data frame to keep only features with adjusted p-value < 0.05
significant_results_dfHC_EPL_numberic  <- subset(results_df, Adjusted_P_Value < 0.05)
# Extract the names of significant features
significant_featuresHC_EPL_numberic <- significant_results_dfHC_EPL_numberic$Variable


columns_to_keep <- c("Stage", significant_featuresHC_EPL_numberic )


filtered_T_Chemicals_dataHC_EPL_numberic <- T_Chemicals_data_HC_EPL2_numberic[, significant_featuresHC_EPL_numberic, drop = FALSE]

# Select these columns from the transposed T_Chemicals_data
# Make sure "Stage" is a column in T_Chemicals_data
if("Stage" %in% colnames(T_Chemicals_data_HC_EPL2_numberic)) {
  filtered_T_Chemicals_dataHC_EPL_numberic <- T_Chemicals_data_HC_EPL2_numberic[, columns_to_keep, drop = FALSE]
} else {
  stop("Column 'Stage' not found in T_Chemicals_data_HC_EPL2_numberic")
}


# for category variables: The following results indicates no significant difference


# Explicitly perform Chi-square tests for each categorical variable

for (col_index in 319:320) {
  # Convert the column to a factor
  T_Chemicals_data_HC_EPL2[, col_index] <- factor(T_Chemicals_data_HC_EPL2[, col_index])
  
  # Ensure that the comparison column is also a factor if it is categorical
  comparison_column <- factor(T_Chemicals_data_HC_EPL2[, 3])  # Replace '2' with the correct column index if needed
  
  # Perform the Chi-square test
  chi_square_result <- tryCatch({
    chisq.test(table(T_Chemicals_data_HC_EPL2[, col_index], comparison_column))
  }, warning = function(w) {
    w
  }, error = function(e) {
    e
  }, finally = {
    cat("Processed column", col_index, "\n")
  })
  
  # Check if result is a chi-square test object and print p-value
  if (inherits(chi_square_result, "htest")) {
    cat("Variable:", colnames(T_Chemicals_data_HC_EPL2)[col_index], "- Chi-square p-value:", chi_square_result$p.value, "\n")
    # Add the Chi-square result to the results list
    results_list[[length(results_list) + 1]] <- list(
      Variable = colnames(T_Chemicals_data_HC_EPL2)[col_index],
      Chi_Square_Test = chi_square_result
    )
  } else {
    cat("Error or warning occurred in Chi-square test for column", col_index, "\n")
  }
}


##------------------------------------------------
## EE vs EPL


EE_EPL_addition = filtered_T_Chemicals_dataEE_EPL_numberic
EE_EPL_addition <-  EE_EPL_addition%>% mutate(Stage = ifelse(Stage =="EE",1,0))
EE_EPL_addition$Stage<-factor(EE_EPL_addition$Stage,levels = c(0,1),labels = c("EPL","EE"))

model_EE_EPL_addition <- glm(Stage ~ .,EE_EPL_addition , family = binomial)
model.both_EE_EPL_addition <- stepAIC(model_EE_EPL_addition, direction = 'both') # 'forward', 'backward'
# 查看模型结果
summary(model.both_EE_EPL_addition)$coefficients
##-------------------------------------------------------------------------------------------------
## cancer vs epl 

## Kw test with BH correction

library(stats)

T_Chemicals_data_CA_EPL <- read.csv("", check.names = FALSE)

# Create an empty list for storing Kruskal-Wallis and BH correction results
results_list <- list()

# Create a vector to store variable names
variable_names <- colnames(T_Chemicals_data_CA_EPL)[3:317]

for (i in 3:317) {
  # Execute Kruskal-Wallis test
  kw_result <- kruskal.test(T_Chemicals_data_CA_EPL[, i], as.factor(T_Chemicals_data_CA_EPL[, 2]))
  
  # Store Kruskal-Wallis test result in the list along with variable name
  result_entry <- list(Variable = variable_names[i - 2], Kruskal_Test = kw_result)
  results_list[[i - 2]] <- result_entry
}

# Extract p-values from the Kruskal-Wallis test results
p_values <- sapply(results_list, function(result) result$Kruskal_Test$p.value)

# Apply BH correction to the p-values
adjusted_p_values <- p.adjust(p_values, method = "BH")


# Add the adjusted p-values to the results list

for (i in 1:length(results_list)) {
  results_list[[i]]$adjusted_p_value <- adjusted_p_values[i]
}

# Convert the list to a data frame for easier writing to file
results_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(Variable = x$Variable, P_Value = x$Kruskal_Test$p.value, Adjusted_P_Value = as.numeric(x$adjusted_p_value))
}))

# Write the data frame to a CSV file
write.csv(results_df, "kw_with_BH_diagnosis_CA_EPL.csv", row.names = FALSE)


# Filter the data frame to keep only features with adjusted p-value < 0.05
significant_results_dfCA_EPL  <- subset(results_df, Adjusted_P_Value < 0.05)
# Extract the names of significant features
significant_featuresCA_EPL <- significant_results_dfCA_EPL$Variable


columns_to_keep <- c("Stage", significant_featuresCA_EPL )


filtered_T_Chemicals_dataCA_EPL <- T_Chemicals_data_CA_EPL[, significant_featuresCA_EPL, drop = FALSE]


# Make sure "Stage" is a column in T_Chemicals_data
if("Stage" %in% colnames(T_Chemicals_data_CA_EPL)) {
  filtered_T_Chemicals_dataCA_EPL <- T_Chemicals_data_CA_EPL[, columns_to_keep, drop = FALSE]
} else {
  stop("Column 'Stage' not found in T_Chemicals_data_CA_EPL")
}

##------------------------------------------------------------
## plot ROC for EPL vs CA

library(stats)

T_Chemicals_data_CA_EPL <- read.csv("", check.names = FALSE)

# Create an empty list for storing Kruskal-Wallis and BH correction results
results_list <- list()

# Create a vector to store variable names
variable_names <- colnames(T_Chemicals_data_CA_EPL)[3:317]

for (i in 3:317) {
  # Execute Kruskal-Wallis test
  kw_result <- kruskal.test(T_Chemicals_data_CA_EPL[, i], as.factor(T_Chemicals_data_CA_EPL[, 2]))
  
  # Store Kruskal-Wallis test result in the list along with variable name
  result_entry <- list(Variable = variable_names[i - 2], Kruskal_Test = kw_result)
  results_list[[i - 2]] <- result_entry
}

# Extract p-values from the Kruskal-Wallis test results
p_values <- sapply(results_list, function(result) result$Kruskal_Test$p.value)

# Apply BH correction to the p-values
adjusted_p_values <- p.adjust(p_values, method = "BH")


# Add the adjusted p-values to the results list

for (i in 1:length(results_list)) {
  results_list[[i]]$adjusted_p_value <- adjusted_p_values[i]
}

# Convert the list to a data frame for easier writing to file
results_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(Variable = x$Variable, P_Value = x$Kruskal_Test$p.value, Adjusted_P_Value = as.numeric(x$adjusted_p_value))
}))

# Write the data frame to a CSV file
write.csv(results_df, "kw_with_BH_diagnosis_CA_EPL.csv", row.names = FALSE)


# Filter the data frame to keep only features with adjusted p-value < 0.05
significant_results_dfCA_EPL  <- subset(results_df, Adjusted_P_Value < 0.05)
# Extract the names of significant features
significant_featuresCA_EPL <- significant_results_dfCA_EPL$Variable


columns_to_keep <- c("Stage", significant_featuresCA_EPL )


filtered_T_Chemicals_dataCA_EPL <- T_Chemicals_data_CA_EPL[, significant_featuresCA_EPL, drop = FALSE]


if("Stage" %in% colnames(T_Chemicals_data_CA_EPL)) {
  filtered_T_Chemicals_dataCA_EPL <- T_Chemicals_data_CA_EPL[, columns_to_keep, drop = FALSE]
} else {
  stop("Column 'Stage' not found in T_Chemicals_data_CA_EPL")
}



library(dplyr) 
library(data.table) 
library(caTools) 
library(pROC) 
library(ggplot2) 
library(ggpubr) 
library(ggprism) 
library(MASS)
library(caret) 
library(boot)



step_CA_EPL <-  step_CA_EPL %>% mutate(Stage = ifelse(Stage =="EPL",1,0))
step_CA_EPL$Stage<-factor(step_CA_EPL$Stage,levels = c(0,1),labels = c("CA","EPL"))
# Load the caret package
library(caret)

set.seed(123)  
split_ratio <- 0.7  
split_vector <- sample.split(step_CA_EPL$Stage, SplitRatio = split_ratio)

train_data <- step_CA_EPL[split_vector, ]  
test_data <- step_CA_EPL[!split_vector, ]  

model <- glm(Stage ~ ., data = train_data, family = binomial)

summary(model)


predictions <- predict(model, newdata = test_data, Stage = "response")


threshold <- 0.5  
predicted_classes <- ifelse(predictions >= threshold, 1, 0)


confusion_matrix <- table(test_data$Stage, predicted_classes)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)


print(confusion_matrix)
cat("Accuracy: ", accuracy, "\n")
cat("Precision: ", precision, "\n")
cat("Recall: ", recall, "\n")
cat("F1 Score: ", f1_score, "\n")

roc_obj <- roc(test_data$Stage, predictions,percent=TRUE, plot=TRUE, ci=TRUE)
roc_obj
roc_auc <- auc(roc_obj)


roc_data <- data.frame(roc_obj$specificities, roc_obj$sensitivities)





p1=ggplot(roc_data, aes(x =  roc_obj.specificities, y = roc_obj.sensitivities)) +
  geom_line(color = "#8092C4", size = 2) +
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), linetype = "dashed", color = "#F47F72",size=1.2) +
  geom_text(aes(x = 0.2, y = 0.2, label = paste("AUC =", round(roc_auc, 3))), size = 6, color = "black") +
  coord_cartesian(xlim = c(1.02, -0.01), ylim = c(-0.01, 1.02),expand=F) +
  theme_minimal()+
  labs(x = "
       Specificity", y = "Sensitivity") +
  ggtitle("EPL vs CA") +
  theme(plot.title = element_text(size = 14, face = "bold"))+
  theme_prism(border = T)



ggsave(filename = "EPL vs CA.pdf",p1, width = 6, height = 5)

##--------------------------------------------------------------------------------------
## plot ROC for EPL vs HC

library(stats)


T_Chemicals_data_HC_EPL <- read.csv("", check.names = FALSE)
results_list2 <- list()

# Create a vector to store variable names
variable_names2 <- colnames(T_Chemicals_data_HC_EPL)[3:317]

for (i in 3:317) {
  # Execute Kruskal-Wallis test
  kw_result <- kruskal.test(T_Chemicals_data_HC_EPL[, i], as.factor(T_Chemicals_data_HC_EPL[, 2]))
  
  # Store Kruskal-Wallis test result in the list along with variable name
  result_entry <- list(Variable = variable_names2[i - 2], Kruskal_Test = kw_result)
  results_list2[[i - 2]] <- result_entry
}

# Extract p-values from the Kruskal-Wallis test results
p_values2 <- sapply(results_list2, function(result) result$Kruskal_Test$p.value)

# Apply BH correction to the p-values
adjusted_p_values <- p.adjust(p_values2, method = "BH")

# Add the adjusted p-values to the results list

for (i in 1:length(results_list2)) {
  results_list2[[i]]$adjusted_p_value <- adjusted_p_values[i]
}

# Convert the list to a data frame for easier writing to file
results_df2 <- do.call(rbind, lapply(results_list2, function(x) {
  data.frame(Variable = x$Variable, P_Value = x$Kruskal_Test$p.value, Adjusted_P_Value = as.numeric(x$adjusted_p_value))
}))

# Write the data frame to a CSV file
write.csv(results_df2, "kw_with_BH_diagnosisHC_EPL.csv", row.names = FALSE)


# Filter the data frame to keep only features with adjusted p-value < 0.05
significant_results_dfHC_EPL  <- subset(results_df2, Adjusted_P_Value < 0.05)
# Extract the names of significant features
significant_featuresHC_EPL <- significant_results_dfHC_EPL$Variable


columns_to_keep <- c("Stage", significant_featuresHC_EPL )


filtered_T_Chemicals_dataHC_EPL <- T_Chemicals_data_HC_EPL[, significant_featuresHC_EPL, drop = FALSE]


if("Stage" %in% colnames(T_Chemicals_data_HC_EPL)) {
  filtered_T_Chemicals_dataHC_EPL <- T_Chemicals_data_HC_EPL[, columns_to_keep, drop = FALSE]
} else {
  stop("Column 'Stage' not found in T_Chemicals_data_HC_EPL")
}


library(dplyr) 
library(data.table) 
library(pROC) 
library(ggplot2) 
library(ggprism) 
library(MASS)
library(caret) 
library(caTools)
library(boot)
library(pROC)
library(ggpubr)

step_HC_EPL <-  step_HC_EPL %>% mutate(Stage = ifelse(Stage =="EPL",1,0))
step_HC_EPL$Stage<-factor(step_HC_EPL$Stage,levels = c(0,1),labels = c("HC","EPL"))
# Load the caret package
library(caret)

set.seed(123)  
split_ratio <- 0.7  
split_vector2 <- sample.split(step_HC_EPL$Stage, SplitRatio = split_ratio)

train_data2 <- step_HC_EPL[split_vector2, ]  
test_data2 <- step_HC_EPL[!split_vector2, ]  


model2 <- glm(Stage ~ ., data = train_data2, family = binomial)

summary(model2)


predictions2 <- predict(model2, newdata = test_data2, Stage = "response")


threshold <- 0.5  
predicted_classes2 <- ifelse(predictions2 >= threshold, 1, 0)


confusion_matrix2 <- table(test_data2$Stage, predicted_classes2)
accuracy2 <- sum(diag(confusion_matrix2)) / sum(confusion_matrix2)
precision2<- confusion_matrix2[2, 2] / sum(confusion_matrix2[, 2])
recall2 <- confusion_matrix2[2, 2] / sum(confusion_matrix2[2, ])
f1_score2 <- 2 * precision2 * recall2 / (precision2 + recall2)


print(confusion_matrix2)
cat("Accuracy: ", accuracy2, "\n")
cat("Precision: ", precision2, "\n")
cat("Recall: ", recall2, "\n")
cat("F1 Score: ", f1_score2, "\n")

roc_obj2 <- roc(test_data2$Stage, predictions2,percent=TRUE, plot=TRUE, ci=TRUE)
roc_obj2
roc_auc2 <- auc(roc_obj2)


roc_data2 <- data.frame(roc_obj2$specificities, roc_obj2$sensitivities)





p2=ggplot(roc_data2, aes(x =  roc_obj2.specificities, y = roc_obj2.sensitivities)) +
  geom_line(color = "#8092C4", size = 2) +
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), linetype = "dashed", color = "#F47F72",size=1.2) +
  geom_text(aes(x = 0.2, y = 0.2, label = paste("AUC =", round(roc_auc2, 3))), size = 6, color = "black") +
  coord_cartesian(xlim = c(1.02, -0.01), ylim = c(-0.01, 1.02),expand=F) +
  theme_minimal()+
  labs(x = "
       Specificity", y = "Sensitivity") +
  ggtitle("EPL vs HC") +
  theme(plot.title = element_text(size = 14, face = "bold"))+
  theme_prism(border = T)

ggsave(filename = "EPL vs HC.pdf",p2, width = 6, height = 5)

## combine plot

p = ggplot() + 
  geom_ribbon(data = roc_data, aes(x = roc_obj.specificities, ymin = 0, ymax = roc_obj.sensitivities), fill = "#E4E7A2", alpha = 0.3) +
  geom_line(data = roc_data, aes(x =  roc_obj.specificities, y = roc_obj.sensitivities),color = "#E4E7A2", size = 1.2) +
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), linetype = "dashed", color = "#F47F72",size=1.2) +
  geom_text(aes(x = 0.2, y = 0.2, label = paste("AUC =", round(roc_auc, 3))), size = 7, color = "#E8BE74") +
  coord_cartesian(xlim = c(1.02, -0.01), ylim = c(-0.01, 1.02),expand=F) +
  theme_minimal()+
  labs(x = "
       Specificity", y = "Sensitivity") +
  ggtitle("ROC curve") +
  theme(plot.title = element_text(size = 14, face = "bold"))+
  theme_prism(border = T)

p = p + 
  geom_ribbon(data = roc_data2, aes(x = roc_obj2.specificities, ymin = 0, ymax = roc_obj2.sensitivities), fill = "#B5CBE2", alpha = 0.3) +
  geom_line(data = roc_data2, aes(x =  roc_obj2.specificities, y = roc_obj2.sensitivities),color = "#B5CBE2", size = 1.2) +
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), linetype = "dashed", color = "#F47F72",size=1.2) +
  geom_text(aes(x = 0.2, y = 0.3, label = paste("AUC =", round(roc_auc2, 3))), size = 7, color = "#55759E") +
  coord_cartesian(xlim = c(1.02, -0.01), ylim = c(-0.01, 1.02),expand=F) +
  theme_minimal()+
  labs(x = "
       Specificity", y = "Sensitivity") +
  ggtitle("ROC curve") +
  theme(plot.title = element_text(size = 14, face = "bold"))+
  theme_prism(border = T)

p

ggsave(filename = "Combine_ROC.pdf",p, width = 6, height = 5)

##-----------------------------------------------------------------------------
#plot circle heatmap for EPL VS CA

numerical_columns <- step_CA_EPL[, sapply(step_CA_EPL, is.numeric)]
scaled_numerical_columns <- scale(numerical_columns)

# Replace the original numerical columns with the scaled ones
step_CA_EPL[, sapply(step_CA_EPL, is.numeric)] <- scaled_numerical_columns

library('ComplexHeatmap')
library('circlize')


madt2<-as.matrix(step_CA_EPL[,-1])

Heatmap(madt2)## plot a normal heatmap

ann_row = data.frame(pathway=c(rep("EPL",14),rep("CA",101)))
row.names(ann_row) = rownames(madt2)
ann_row <- as.matrix(ann_row)


range(madt2)
mycol=colorRamp2(c(-1.7, 0.3, 2.3),c("#57ab81", "white", "#ff9600"))
circos.heatmap(madt2,split=ann_row,col=mycol)
circos.clear()

#to make the heatmap look good
circos.par(gap.after=c(10))
circos.heatmap(madt2,col=mycol,split=ann_row,rownames.side="outside",
               rownames.col="black",
               rownames.cex=0.7,
               rownames.font=0.1)
circos.clear()



circos.par(gap.after=c(10))
circos.heatmap(madt2,split=ann_row,show.sector.labels = T,col=mycol,rownames.side="outside",track.height = 0.5,
               rownames.col="black",
               rownames.cex=0.7,
               rownames.font=0.1)



## add column name
circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
  if(CELL_META$sector.numeric.index==1){
    cn=colnames(madt2)
    n=length(cn)
    circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(0.8,"mm"),#x坐标
                (1:n)*1.4,
                cn,cex=1,adj=c(0,0.4),facing="inside")
  }
},bg.border=NA)


lgd=Legend(title="z score",col_fun=mycol,direction = c("vertical"))

draw(lgd, x= unit(1.1, "snpc"),just ="rightbottom")

circos.clear()

##--------------------------------------------------------------------------
#plot circle heatmap for EPL VS HC

numerical_columns <- step_HC_EPL[, sapply(step_HC_EPL, is.numeric)]
scaled_numerical_columns <- scale(numerical_columns)

# Replace the original numerical columns with the scaled ones
step_HC_EPL[, sapply(step_HC_EPL, is.numeric)] <- scaled_numerical_columns

library('ComplexHeatmap')
library('circlize')


madt<-as.matrix(step_HC_EPL[,-1])

Heatmap(madt)## plot a normal heatmap

ann_row = data.frame(pathway=c(rep("EPL",14),rep("HC",43)))
row.names(ann_row) = rownames(madt)
ann_row <- as.matrix(ann_row)


range(madt)
mycol=colorRamp2(c(-1.7, 0.3, 2.3),c("#57ab81", "white", "#ff9600"))
circos.heatmap(madt,split=ann_row,col=mycol)
circos.clear()

#to make the heatmap look good
circos.par(gap.after=c(10))
circos.heatmap(madt,col=mycol,split=ann_row,rownames.side="outside",
               rownames.col="black",
               rownames.cex=0.7,
               rownames.font=0.1)
circos.clear()


circos.par(gap.after=c(10))
circos.heatmap(madt,split=ann_row,show.sector.labels = T,col=mycol,rownames.side="outside",track.height = 0.5,
               rownames.col="black",
               rownames.cex=0.7,
               rownames.font=0.1)


## add column name
circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
  if(CELL_META$sector.numeric.index==1){
    cn=colnames(madt)
    n=length(cn)
    circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(0.8,"mm"),
                (1:n)*1.2,
                cn,cex=1,adj=c(0,0.4),facing="inside")
  }
},bg.border=NA)



lgd=Legend(title="z score",col_fun=mycol,direction = c("vertical"))

draw(lgd, x= unit(1.1, "snpc"),just ="rightbottom")

circos.clear()

##-------------------------------------------------------------------
## HC VS ESCC

## find biomarkers bewtween HC and ESCC

## 1. differential chemcais

library(stats)

T_Chemicals_data_HC_ESCC <- read.csv("", check.names = FALSE)

# Create an empty list for storing Kruskal-Wallis and BH correction results
results_list <- list()

# Create a vector to store variable names
variable_names <- colnames(T_Chemicals_data_HC_ESCC)[3:317]

for (i in 3:317) {
  # Execute Kruskal-Wallis test
  kw_result <- kruskal.test(T_Chemicals_data_HC_ESCC[, i], as.factor(T_Chemicals_data_HC_ESCC[, 2]))
  
  # Store Kruskal-Wallis test result in the list along with variable name
  result_entry <- list(Variable = variable_names[i - 2], Kruskal_Test = kw_result)
  results_list[[i - 2]] <- result_entry
}

# Extract p-values from the Kruskal-Wallis test results
p_values <- sapply(results_list, function(result) result$Kruskal_Test$p.value)

# Apply BH correction to the p-values
adjusted_p_values <- p.adjust(p_values, method = "BH")


# Add the adjusted p-values to the results list

for (i in 1:length(results_list)) {
  results_list[[i]]$adjusted_p_value <- adjusted_p_values[i]
}

# Convert the list to a data frame for easier writing to file
results_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(Variable = x$Variable, P_Value = x$Kruskal_Test$p.value, Adjusted_P_Value = as.numeric(x$adjusted_p_value))
}))

# Write the data frame to a CSV file
write.csv(results_df, "kw_with_BH_diagnosis_HC_ESCC.csv", row.names = FALSE)


# Filter the data frame to keep only features with adjusted p-value < 0.05
significant_results_dfHC_ESCC  <- subset(results_df, Adjusted_P_Value < 0.05)
# Extract the names of significant features
significant_featuresHC_ESCC <- significant_results_dfHC_ESCC$Variable


columns_to_keep <- c("Stage", significant_featuresHC_ESCC )


filtered_T_Chemicals_dataHC_ESCC <- T_Chemicals_data_HC_ESCC[, significant_featuresHC_ESCC, drop = FALSE]

# Select these columns from the transposed T_Chemicals_data
# Make sure "Stage" is a column in T_Chemicals_data
if("Stage" %in% colnames(T_Chemicals_data_HC_ESCC)) {
  filtered_T_Chemicals_dataHC_ESCC <- T_Chemicals_data_HC_ESCC[, columns_to_keep, drop = FALSE]
} else {
  stop("Column 'Stage' not found in T_Chemicals_data_HC_ESCC")
}


## 2. stepwise regression selection 


library(dplyr) 
library(data.table)
library(caTools) 
library(pROC) 
library(ggplot2) 
library(ggpubr) 
library(ggprism) 
library(MASS)
library(caret) 
library(boot)


filtered_T_Chemicals_dataHC_ESCC <- filtered_T_Chemicals_dataHC_ESCC%>% mutate(Stage = ifelse(Stage =="EC",1,0))
filtered_T_Chemicals_dataHC_ESCC$Stage<-factor(filtered_T_Chemicals_dataHC_ESCC$Stage,levels = c(0,1),labels = c("hlthy","EC"))

model_HC_ESCC <- glm(Stage ~ .,filtered_T_Chemicals_dataHC_ESCC, family = binomial)
model.both_HC_ESCC <- stepAIC(model_HC_ESCC, direction = 'both') # 'forward', 'backward'

summary(model.both_HC_ESCC)$coefficients


## 3. model performance


step_HC_ESCC <-  step_HC_ESCC %>% mutate(Stage = ifelse(Stage =="EC",1,0))
step_HC_ESCC$Stage<-factor(step_HC_ESCC$Stage,levels = c(0,1),labels = c("hlthy","EC"))
# Load the caret package
library(caret)

set.seed(123)  
split_ratio <- 0.7 
split_vector <- sample.split(step_HC_ESCC$Stage, SplitRatio = split_ratio)

train_data <- step_HC_ESCC[split_vector, ]  
test_data <- step_HC_ESCC[!split_vector, ]  


model <- glm(Stage ~ ., data = train_data, family = binomial)

summary(model)



predictions <- predict(model, newdata = test_data,type = "response")


threshold <- 0.5  
predicted_classes <- ifelse(predictions >= threshold, 1, 0)


confusion_matrix <- table(test_data$Stage, predicted_classes)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)


print(confusion_matrix)
cat("Accuracy: ", accuracy, "\n")
cat("Precision: ", precision, "\n")
cat("Recall: ", recall, "\n")
cat("F1 Score: ", f1_score, "\n")


roc_obj <- roc(test_data$Stage, predictions,percent=TRUE, plot=TRUE, ci=TRUE)
roc_obj
roc_auc <- auc(roc_obj)


roc_data <- data.frame(roc_obj$specificities, roc_obj$sensitivities)


## #plot circle heatmap for ESCC VS HC

step_HC_ESCC<-filtered_T_Chemicals_dataHC_ESCC[,c(1,3,121)]


numerical_columns <- step_HC_ESCC[, sapply(step_HC_ESCC, is.numeric)]
scaled_numerical_columns <- scale(numerical_columns)

# Replace the original numerical columns with the scaled ones
step_HC_ESCC[, sapply(step_HC_ESCC, is.numeric)] <- scaled_numerical_columns

library('ComplexHeatmap')
library('circlize')


madt<-as.matrix(step_HC_ESCC[,-1])

Heatmap(madt)## plot a normal heatmap

ann_row = data.frame(pathway=c(rep("ESCC",101),rep("HC",43)))
row.names(ann_row) = rownames(madt)
ann_row <- as.matrix(ann_row)


range(madt)
mycol=colorRamp2(c(-1.7, 0.3, 2.3),c("#57ab81", "white", "#ff9600"))
circos.heatmap(madt,split=ann_row,col=mycol)
circos.clear()

#to make the heatmap look good
circos.par(gap.after=c(10))
circos.heatmap(madt,col=mycol,split=ann_row,rownames.side="outside",
               rownames.col="black",
               rownames.cex=0.7,
               rownames.font=0.1)
circos.clear()



circos.par(gap.after=c(10))
circos.heatmap(madt,split=ann_row,show.sector.labels = T,col=mycol,rownames.side="outside",track.height = 0.5,
               rownames.col="black",
               rownames.cex=0.7,
               rownames.font=0.1)


## add column name
circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
  if(CELL_META$sector.numeric.index==1){
    cn=colnames(madt)
    n=length(cn)
    circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(0.8,"mm"),#x坐标
                (1:n)*1.2,#y坐标
                cn,cex=1,adj=c(0,0.4),facing="inside")
  }
},bg.border=NA)



lgd=Legend(title="z score",col_fun=mycol,direction = c("vertical"))

draw(lgd, x= unit(1.1, "snpc"),just ="rightbottom")

circos.clear()

