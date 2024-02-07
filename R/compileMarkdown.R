#' compileMarkdown
#'
#' @description This function compiles the markdown file with the analysis results
#'
#' @param lassoPredictorsPath String | Path to the mean number of predictors selected by Lasso in each generation.
#' @param geneticAlgorithm Genetic algorithm object.
#' @param originalDiagnosis Original diagnostics of the patients.
#' @param clinicData Dataset of clinic data that will be used.
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::compileMarkdown(savingName = savingName, mlAlgorithm = mlAlgorithm, lassoPredictorsPath = lassoPredictorsPath, geneticAlgorithm = geneticAlgorithm, clinicData = clinicData, classVariable = classVariable)
#'

compileMarkdown <- function(savingName, mlAlgorithm, lassoPredictorsPath, geneticAlgorithm, originalDiagnosis, clinicData, classVariable){

  name <- paste("GA", savingName, sep="_")

  dirPath <- paste(savingName, "analysisData", name, sep = "/")

  pcaPath <- paste(dirPath, "PCA.tsv", sep="_")
  pca <- read.table(pcaPath, header = TRUE, sep = "\t", row.names = 1)

  gaPathNumeric <- paste(dirPath, "NumericTable.tsv", sep="_")
  gaPathTotal <- paste(dirPath, "TotalTable.tsv", sep="_")

  numeric <- read.table(gaPathNumeric, header = TRUE, sep = "\t", row.names = 1)
  total <- read.table(gaPathTotal, header = TRUE, sep = "\t", row.names = 1)

  algorithmPrecisions <- geneticAlgorithm@summary[,1]

  outputName <- paste("analysisResult_", name, ".html", sep = "")

  outputPath <- paste("./", savingName, "/", sep = "")

  lassoPredictors <- NULL

  if(mlAlgorithm == "Lasso"){

    if(justAnalysis){
      lassoPredictors <- readRDS(lassoPredictorsPath)
    } else {

      gaPath <- paste(dirPath, "Predictors.rds", sep="_")
      lassoPredictors <- readRDS(gaPath)
    }


  }

  rmarkdown::render(input = system.file("data", "analysisResult.Rmd", package = "MLASDO"),
                   params = list(
                                 algorithmPrecisions = algorithmPrecisions,
                                 lassoPredictors = lassoPredictors,
                                 originalDiagnosis = originalDiagnosis,
                                 mlAlgorithm = mlAlgorithm,
                                 clinicData = clinicData,
                                 pcaAnalisis = pca,
                                 numericTable = numeric,
                                 totalTable = total,
                                 classVariable = classVariable),
                   output_file = outputName,
                   output_dir = outputPath
                   )
}
