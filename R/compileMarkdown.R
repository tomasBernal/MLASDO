#' compileMarkdown
#'
#' This function compiles the markdown file with the analysis results
#'
#' @param savingName Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#' @param clinicData Dataset of clinic data that will be used.
#' @param classVariable Target variable, which must be binary, meaning it has two possible values.
#'
#' @export
#'
#' @examples
#'
#' MLASDO::compileMarkdown(savingName = savingName, clinicData = clinicData, classVariable = classVariable)
#'

compileMarkdown <- function(savingName, clinicData, classVariable){

  name <- paste("GA", savingName, sep="_")

  pcaPath <- paste(name, "PCA.tsv", sep="_")
  pca <- read.table(pcaPath, header = TRUE, sep = "\t", row.names = 1)

  gaPathNumeric <- paste(name, "NumericTable.tsv", sep="_")
  gaPathTotal <- paste(name, "TotalTable.tsv", sep="_")

  numeric <- read.table(gaPathNumeric, header = TRUE, sep = "\t", row.names = 1)
  total <- read.table(gaPathTotal, header = TRUE, sep = "\t", row.names = 1)

  rmarkdown::render(input = system.file("data", "analysisResult.Rmd", package = "MLASDO"),
                   params = list(
                                 clinicData = clinicData,
                                 pcaAnalisis = pca,
                                 numericTable = numeric,
                                 totalTable = total,
                                 classVariable = classVariable),
                   output_file = "analysisResult.html",
                   output_dir = "./")
}
