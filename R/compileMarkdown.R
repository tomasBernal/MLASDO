#' compileMarkdown
#'
#' @description This function compiles the markdown file with the analysis results
#'
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
#' MLASDO::compileMarkdown(savingName = savingName, geneticAlgorithm = geneticAlgorithm, clinicData = clinicData, classVariable = classVariable)
#'

compileMarkdown <- function(savingName, geneticAlgorithm, originalDiagnosis, clinicData, classVariable){

  name <- paste("GA", savingName, sep="_")

  pcaPath <- paste(name, "PCA.tsv", sep="_")
  pca <- read.table(pcaPath, header = TRUE, sep = "\t", row.names = 1)

  gaPathNumeric <- paste(name, "NumericTable.tsv", sep="_")
  gaPathTotal <- paste(name, "TotalTable.tsv", sep="_")

  numeric <- read.table(gaPathNumeric, header = TRUE, sep = "\t", row.names = 1)
  total <- read.table(gaPathTotal, header = TRUE, sep = "\t", row.names = 1)

  algorithmPrecisions <- ga@summary[,1]

  rmarkdown::render(input = system.file("data", "analysisResult.Rmd", package = "MLASDO"),
                   params = list(
                                 algorithmPrecisions = algorithmPrecisions,
                                 originalDiagnosis = originalDiagnosis,
                                 clinicData = clinicData,
                                 pcaAnalisis = pca,
                                 numericTable = numeric,
                                 totalTable = total,
                                 classVariable = classVariable),
                   output_file = "analysisResult.html",
                   output_dir = "./")
}
