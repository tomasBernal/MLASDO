#' performPCAAnalysis
#'
#' This function performs PCA analysis on the changed diagnoses after the execution of the genetic algorithm.
#'
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#' @param mlAlgorithm String | Machine Learning algorithm to be applied, the options are: Lasso or RF (Random Forest).
#'
#' @param changedOmicData Dataset | Dataset of omic data that will be used.
#' @param selectedData Dataset | Dataset of omic data with only the predictors selected by the Lasso model.
#'
#' @param idColumn String | Variable that indicates the identifier of each patient in both datasets. If the user does not specify a path to his own data, the value for the sample data, Trial, will be used.
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::performPCAAnalysis(savingName = savingName, mlAlgorithm = mlAlgorithm, changedOmicData = changedOmicData, selectedData = selectedData, idColumn = idColumn, classVariable = classVariable)


performPCAAnalysis <- function(
    savingName,
    mlAlgorithm,
    changedOmicData,
    selectedData,
    idColumn,
    classVariable
){

  #### REQUIRED LIBRARIES ####
  library(dplyr) # For select function

  #### DATA READING ####
  omic <- changedOmicData

  omicShow <- changedOmicData

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

  vExp <- omicResult$sdev^2

  pExp <- vExp / sum(vExp)
  pExp <- round(pExp, 2)


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


  selectedData[[idColumn]] <- NULL

  # We need to remove non numeric columns
  selectedData <- dplyr::select(selectedData, where(is.numeric))

  # We need to remove constant or zero variance columns
  uselessColumns <- sapply(selectedData, function(x) length(unique(x)) == 1)

  selectedData <- selectedData[, !uselessColumns]

  # Performing prcomp which gives us the deviations of the principal components
  omicResult <- prcomp(selectedData, scale. = TRUE)

  vExpSelected <- omicResult$sdev^2

  pExpSelected <- vExpSelected / sum(vExpSelected)
  pExpSelected <- round(pExpSelected, 2)

  # I convert the PCA analysis into a data frame and add the column indicating the diagnosis.
  pca <- as.data.frame(omicResult$x, stringsAsFactors=F)

  # I add the changed diagnoses to the PCA data frame.
  pca <- cbind(classVariable = omicShow[[classVariable]], pca)

  # Changing the column name
  colnames(pca)[colnames(pca) == "classVariable"] <- classVariable

  # Save the pca analysis
  gaPath <- paste(dirPath, "PCA_Selected.tsv", sep="_")
  write.table(pca, gaPath, row.names = T, col.names = T, sep =  '\t')


  pcaInfo <- data.frame(
    PC1 = pExp[1],
    PC2 = pExp[2],
    PC1Selected = pExpSelected[1],
    PC2Selected = pExpSelected[2]
  )

  # Save the pca analysis
  gaPath <- paste(dirPath, "PCA_Info.tsv", sep="_")
  write.table(pcaInfo, gaPath, row.names = T, col.names = T, sep =  '\t')
}
