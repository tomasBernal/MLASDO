#' performPCAAnalysis
#'
#' This function performs PCA analysis on the changed diagnoses after the execution of the genetic algorithm.
#'
#'
#' @param idColumn Variable that indicates the identifier of each patient in both datasets.
#' @param omicData Dataset of omic data that will be used.
#' @param activePredictors Predictors on which the study of the ratios will be conducted after the genetic algorithm has been performed.
#' @param classVariable Target variable, which must be binary, meaning it has two possible values.
#' @param savingName Name under which the model and solution were saved after the execution of the genetic algorithm.
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::performPCAAnalysis(idColumn = idColumn, omicData = omicData, activePredictors = activePredictors, classVariable = classVariable, savingName = savingName)


performPCAAnalysis <- function(
    idColumn,
    omicData,
    activePredictors,
    classVariable,
    savingName
){

  #### REQUIRED LIBRARIES ####
  library(dplyr) # For select function

  #### DATA READING ####
  omic <- omicData

  #### GENETIC ALGORITHM SOLUTION READING ####
  name <- paste("GA", savingName, sep="_")
  solutionPath <- paste(name, "Solution.rds", sep="_")
  solutionGA <- readRDS(solutionPath)

  #### PCA ANALYSIS ####

  # We need to remove the column that indicates the patient id
  omic[[idColumn]] <- NULL

  # We need to remove non numeric columns
  omic <- dplyr::select(omic, where(is.numeric))

  # We need to remove constant or zero variance columns
  uselessColumns <- sapply(omic, function(x) length(unique(x)) == 1)

  omic <- omic[, !uselessColumns]

  # Performing prcomp which gives us the deviations of the principal components
  omicResult <- prcomp(omic, scale. = TRUE)

  # Creating a copy of the original diagnosis
  changedDiagnoses <- omicData[[classVariable]]

  # I create two text strings to represent the changed diagnoses
  firstGroup = paste("Case", "Control", sep = "2")
  secondGroup = paste("Control", "Case", sep = "2")

  # Converting original diagnoses from strings to integers
  changedDiagnoses <- ifelse(changedDiagnoses == "Case", 1, 0)

  # I apply the solution from the genetic algorithm to the original diagnoses
  # I use an XOR function for this purpose
  changedDiagnoses <- bitwXor(changedDiagnoses, solutionGA)

  # I obtain the indices where the genetic algorithm solution indicates changes in someone's diagnosis
  changeIndices <- which(solutionGA == 1)

  # Iterating through the indices that have been modified
  for(i in changeIndices){

    # If the gene, after modification, had the value 1
    # It means it had the value 0 (Control)
    if(changedDiagnoses[i] == 1){

      # Therefore, I set it to 2
      changedDiagnoses[i] <- secondGroup

      # If not, it means it had the value 1 (Case)
    } else{

      # Therefore, I set it to 3
      changedDiagnoses[i] <- firstGroup
    }
  }

  changedDiagnoses <- ifelse(changedDiagnoses == 0, "Control", changedDiagnoses)
  changedDiagnoses <- ifelse(changedDiagnoses == 1, "Case", changedDiagnoses)

  # After this, I have the following values:
  #     0 -> Gene remained with the value 0 of the variable being predicted
  #     1 -> Gene remained with the value 1 of the variable being predicted
  #     2 -> Gene changed from the value 0 to the value 1
  #     3 -> Gene changed from the value 1 to the value 0


  # I convert the PCA analysis into a data frame and add the column indicating the diagnosis.
  pca <- as.data.frame(omicResult$x, stringsAsFactors=F)

  # I add the changed diagnoses to the PCA data frame.
  pca <- cbind(classVariable = changedDiagnoses, pca)

  # Changing the column name
  colnames(pca)[colnames(pca) == "classVariable"] <- classVariable

  # Save the pca analysis
  solutionPath <- paste(name, "PCA.tsv", sep="_")
  write.table(pca, solutionPath, row.names = T, col.names = T, sep =  '\t')
}
