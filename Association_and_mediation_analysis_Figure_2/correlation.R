# Load necessary libraries
library(tidyverse)
library(nnet)

NEM<- read.csv("", check.names = T)
EM <- read.csv("", check.names = T)
org_endo <- read.csv("", check.names = T)
org_exo <- read.csv("", check.names = T)
habit <- read.csv("", check.names = T)
health_outcome <- read.csv("", check.names = T)
covariates <- read.csv("", check.names = T)

# Convert specific columns to factors
covariates[, c(2,3,6,7)] <- lapply(covariates[, c(2,3,6,7)], as.factor)
habit[, 2:10] <- lapply(habit[, 2:10], as.factor)


# Perform z-score transformation on  columns

EM[, 2:10] <- scale(EM[, 2:10])
NEM[, 2:8] <- scale(NEM[, 2:8])
org_endo[, 2:228] <- scale(org_endo[, 2:228])
org_exo[, 2:72] <- scale(org_exo[, 2:72])


NEM<- NEM[,-1]
EM <- EM[,-1]
org_endo <- org_endo[,-1]
org_exo <-org_exo[,-1]
habit <- habit[,-1]
health_outcome <- health_outcome[,-1]
covariates <- covariates[,-1]



##---------------------------------------------------------------------------------------------------------------------------------
## define function

## define function

perform_analysis <- function(chemical_data, habit_data, health_data, covariate_data) {
  results <- list()
  
  # Correlation with habits (multinomial logistic regression)
  
  for (habit in colnames(habit_data)) {
    for (chemical in colnames(chemical_data)) {
      formula <- as.formula(paste(habit, "~", chemical, "+", paste(colnames(covariate_data), collapse="+")))
      # Check if habit has more than two levels
      if (is.factor(habit_data[[habit]]) && nlevels(habit_data[[habit]]) > 2) {
        model <- multinom(formula, data = cbind(habit_data, chemical_data, covariate_data))
      } else {
        model <- glm(formula, data = cbind(habit_data, chemical_data, covariate_data), family = "binomial")
      }
      results[[paste("Habit", habit, "Chemical", chemical, sep = "_")]] <- summary(model)
    }
  }
  
  
  # Correlation with health outcomes (linear regression)
  for (outcome in colnames(health_data)) {
    for (chemical in colnames(chemical_data)) {
      formula <- as.formula(paste(outcome, "~", chemical, "+", paste(colnames(covariate_data), collapse="+")))
      model <- glm(formula, data = cbind(health_data, chemical_data, covariate_data))
      results[[paste("HealthOutcome", outcome, "Chemical", chemical, sep = "_")]] <- summary(model)
    }
  }
  
  return(results)
}


##----------------------------------------------------------------------------------------------------------------
## Run analysis for each chemical data

results_NEMs <- perform_analysis(NEM, habit, health_outcome, covariates)

results_EMs <- perform_analysis(EM, habit, health_outcome, covariates)
results_org_endo <- perform_analysis(org_endo, habit, health_outcome, covariates)
results_org_exo <- perform_analysis(org_exo, habit, health_outcome, covariates)





##--------------------------------------------------------------------------------------------------------------------------
##  output


## extraction function
extract_model_info <- function(model) {
  # Assuming 'model' is a glm or multinom model object
  # Modify this function based on the specifics of what you want to extract
  coef_summary <- summary(model)$coefficients  # Extracting coefficients
  
  # Convert to dataframe and add any other needed information
  df <- as.data.frame(coef_summary)
  return(df)
}


## Applying the Extraction Function to All Results

combine_results <- function(results) {
  all_results_df <- data.frame()
  
  for (name in names(results)) {
    model_info <- extract_model_info(results[[name]])
    model_info$Model <- name  # Add a column to identify the model
    all_results_df <- rbind(all_results_df, model_info)
  }
  
  return(all_results_df)
}


## Combine and Export Results for Each Chemical Category

# Assuming you have results_nems, results_ems, results_org_endo, results_org_exo from your previous analysis
all_nems <- combine_results(results_NEMs)

write.csv(all_nems, "nems_results.csv", row.names = FALSE)

all_ems <- combine_results(results_ems)
write.csv(all_ems, "ems_results.csv", row.names = FALSE)

all_org_endo <- combine_results(results_org_endo)
write.csv(all_org_endo, "org_endo_results.csv", row.names = FALSE)

all_org_exo <- combine_results(results_org_exo)
write.csv(all_org_exo, "org_exo_results.csv", row.names = FALSE)

#--------------------------------------------------------------------------------------------------------------------------------
## plot correlation chord diagram

library(circlize)

mat1<- read.csv("", check.names = T)

mat1<- mat1[,-2]

col_fun <- colorRamp2(c(-1, 0, 1), c("lightblue", "white", "lightpink"), transparency = 0.1)

set.seed(112)

pdf("chord_cor_habits.pdf", width = 9.5,height = 5)

# circos.clear()

circos.par(start.degree = 90)

chordDiagram(mat1,
             big.gap = 20,
             col = col_fun,
             annotationTrack = "grid",
             preAllocateTracks = list(
               track.height = max(strwidth(unlist(dimnames(mat1))))))


##circos.track(track.index = 1, panel.fun = function(x, y) {
  ##circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              ##facing = "clockwise", niceFacing = TRUE, adj = c(0, 1),
              ## cex = 0.8)}, bg.border = NA)##


highlight.sector(mat1$habit, track.index = 1, col = "#fc82b2",
                 text = "Habits", cex = 1, text.col = "white",
                 niceFacing = TRUE,padding = c(-0.01, 0, -0.9, 0))

#highlight.sector(mat1$chemical_category, track.index = 1, col = "#a3c6e7",
                 #text = "Chemicals", cex = 1, text.col = "white",
                 #niceFacing = TRUE, padding = c(-.1, 0, -.1, 0))

library(ComplexHeatmap)

lgd_links = Legend(at = c(-2.5,  0, 2.5), col_fun = col_fun,
                   title_position = "topleft", title = "cor")

draw(lgd_links, x = unit(1, "npc") - unit(30, "mm"), y = unit(90, "mm"),
     just = c("right"))


dev.off()

circos.clear()
##----------------
# plot2 

library(circlize)

mat2<- read.csv("", check.names = T)


col_fun <- colorRamp2(c(-1, 0, 1), c("lightblue", "white", "lightpink"), transparency = 0.1)

set.seed(112)

pdf("chord_cor_health.pdf", width = 9.5,height = 5)

# circos.clear()

circos.par(start.degree = 90)

chordDiagram(mat2,
             big.gap = 20,
             col = col_fun,
             annotationTrack = "grid",
             preAllocateTracks = list(
               track.height = max(strwidth(unlist(dimnames(mat2))))))


#circos.track(track.index = 1, panel.fun = function(x, y) {
  #circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
             # facing = "clockwise", niceFacing = TRUE, adj = c(0, 1),
             # cex = 0.8 )}, bg.border = NA)


highlight.sector(mat2$health, track.index = 1, col = "#fc82b2",
                 text = "health status indicators", cex = 1, text.col = "white",
                 niceFacing = TRUE, padding = c(-0.01, 0, -0.9, 0))

#highlight.sector(mat2$chemical, track.index = 1, col = "#a3c6e7",
                 #text = "Chemicals", cex = 1, text.col = "white",
                 #niceFacing = TRUE, padding = c(-.5, 0, -.2, 0))

library(ComplexHeatmap)

lgd_links = Legend(at = c(-2.5,  0, 2), col_fun = col_fun,
                   title_position = "topleft", title = "cor")

draw(lgd_links, x = unit(1, "npc") - unit(30, "mm"), y = unit(90, "mm"),
     just = c("right"))


dev.off()




