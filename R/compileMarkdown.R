#' compileMarkdown
#'
#' @description This function compiles the markdown file with the analysis results
#'
#' @param justAnalysis Bool | Indicates whether to perform the analysis directly (TRUE) or to run the genetic algorithm (FALSE). Default value: FALSE.
#' @param lassoPredictorsPath String | Path to the mean number of predictors selected by Lasso in each generation.
#' @param baselinePrecision Decimal | Baseline precision obtained with the model before the detection.
#' @param baselinePredictors Integer | Baseline predictors selected with the model before the detection.
#' @param geneticAlgorithm Array of Strings | Genetic algorithm object.
#' @param originalDiagnosis Array of Strings | Original diagnostics of the patients.
#' @param clinicData Data | Dataset of clinic data that will be used.
#' @param selectedData Data | Dataset of omic data with only the predictors selected by the Lasso model.
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::compileMarkdown(savingName = savingName, mlAlgorithm = mlAlgorithm, lassoPredictorsPath = lassoPredictorsPath, baselinePrecision = baselinePrecision, baselinePredictors = baselinePredictors, geneticAlgorithm = geneticAlgorithm, clinicData = clinicData, selectedData = selectedData, classVariable = classVariable)
#'

compileMarkdown <- function(
    savingName,
    justAnalysis,
    lassoPredictorsPath,
    baselinePrecision,
    baselinePredictors,
    mlAlgorithm,
    geneticAlgorithm,
    originalDiagnosis,
    clinicData,
    selectedData,
    classVariable
    ){

  name <- paste("GA", savingName, sep="_")

  dirPath <- paste(savingName, "analysisData", name, sep = "/")

  pcaPath <- paste(dirPath, "PCA.tsv", sep="_")
  pca <- read.table(pcaPath, header = TRUE, sep = "\t", row.names = 1)

  pcaLasso <- NULL

  if(mlAlgorithm == "Lasso"){

    pcaPath <- paste(dirPath, "PCA_Lasso.tsv", sep="_")
    pcaLasso <- read.table(pcaPath, header = TRUE, sep = "\t", row.names = 1)
  }

  gaPathNumeric <- paste(dirPath, "NumericTable.tsv", sep="_")
  gaPathTotal <- paste(dirPath, "TotalTable.tsv", sep="_")

  numeric <- read.table(gaPathNumeric, header = TRUE, sep = "\t", row.names = 1)
  total <- read.table(gaPathTotal, header = TRUE, sep = "\t", row.names = 1)

  outputName <- paste("analysisResult_", name, ".html", sep = "")

  outputPath <- paste("./", savingName, "/", sep = "")

  lassoPredictors <- NULL

  if(mlAlgorithm == "Lasso"){

    if(justAnalysis){

      lassoPredictors <- readRDS(lassoPredictorsPath)

    } else {

      dirPath <- paste(savingName, "geneticAlgorithm", name, sep = "/")
      gaPath <- paste(dirPath, "Predictors.rds", sep="_")
      lassoPredictors <- readRDS(gaPath)
    }
  }

  rmarkdown::render(input = system.file("data", "analysisResult.Rmd", package = "MLASDO"),
                   params = list(
                                 geneticAlgorithm = geneticAlgorithm,
                                 lassoPredictors = lassoPredictors,
                                 baselinePrecision = baselinePrecision,
                                 baselinePredictors = baselinePredictors,
                                 originalDiagnosis = originalDiagnosis,
                                 mlAlgorithm = mlAlgorithm,
                                 clinicData = clinicData,
                                 selectedData = selectedData,
                                 pcaAnalisis = pca,
                                 pcaAnalisisLasso = pcaLasso,
                                 numericTable = numeric,
                                 totalTable = total,
                                 classVariable = classVariable
                                 ),
                   output_file = outputName,
                   output_dir = outputPath
                   )
}
