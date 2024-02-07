#' performPCAAnalysis
#'
#' @description This function performs PCA analysis on the changed diagnoses after the execution of the genetic algorithm.
#'
#' @param solution Solution obtained after the detection.
#' @param model Best model obtained after the detection.
#'
#' @param idColumn String | Variable that indicates the identifier of each patient in both datasets. If the user does not specify a path to his own data, the value for the sample data, Trial, will be used.
#' @param omicData Dataset of omic data that will be used.
#' @param subsetTrain Subset of train samples.
#' @param activePredictors Array of Strings | Predictors on which the study of the ratios will be conducted after the genetic algorithm has been performed. Default value: All the predictors in clinic data, except classVariable and idColumn.
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::performPCAAnalysis(solution = solution, model = model, idColumn = idColumn, omicData = omicData, subsetTrain = subsetTrain, activePredictors = activePredictors, classVariable = classVariable, savingName = savingName)


performPCAAnalysis <- function(
    solution,
    model,
    idColumn,
    omicData,
    subsetTrain,
    activePredictors,
    classVariable,
    savingName
){

  #### REQUIRED LIBRARIES ####
  library(dplyr) # For select function

  #### DATA READING ####
  omic <- omicData

  omicShow <- omicData

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
  changedDiagnoses <- omicData[subsetTrain,][[classVariable]]

  # I create two text strings to represent the changed diagnoses
  firstGroup = paste("Case", "Control", sep = "2")
  secondGroup = paste("Control", "Case", sep = "2")

  # Converting original diagnoses from strings to integers
  changedDiagnoses <- ifelse(changedDiagnoses == "Case", 1, 0)

  # I apply the solution from the genetic algorithm to the original diagnoses
  # I use an XOR function for this purpose
  changedDiagnoses <- bitwXor(changedDiagnoses, solution)

  # I obtain the indices where the genetic algorithm solution indicates changes in someone's diagnosis
  changeIndices <- which(solution == 1)

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

  omicShow[subsetTrain,][[classVariable]] <- changedDiagnoses

  # I convert the PCA analysis into a data frame and add the column indicating the diagnosis.
  pca <- as.data.frame(omicResult$x, stringsAsFactors=F)

  # I add the changed diagnoses to the PCA data frame.
  pca <- cbind(classVariable = omicShow[[classVariable]], pca)

  # Changing the column name
  colnames(pca)[colnames(pca) == "classVariable"] <- classVariable

  name <- paste("GA", savingName, sep="_")
  dirPath <- paste(savingName, "analysisData", name, sep = "/")

  # Save the pca analysis
  gaPath <- paste(dirPath, "PCA.tsv", sep="_")
  write.table(pca, gaPath, row.names = T, col.names = T, sep =  '\t')
}
