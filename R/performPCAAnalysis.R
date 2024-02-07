#' performPCAAnalysis
#'
#' @description This function performs PCA analysis on the changed diagnoses after the execution of the genetic algorithm.
#'
#' @param model ML Model | Best model obtained after the detection.
#' @param idColumn String | Variable that indicates the identifier of each patient in both datasets. If the user does not specify a path to his own data, the value for the sample data, Trial, will be used.
#' @param changedOmicData Data | Dataset of omic data that will be used.
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::performPCAAnalysis(model = model, idColumn = idColumn, changedOmicData = changedOmicData, classVariable = classVariable, savingName = savingName)


performPCAAnalysis <- function(
    model,
    idColumn,
    changedOmicData,
    classVariable,
    savingName
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
