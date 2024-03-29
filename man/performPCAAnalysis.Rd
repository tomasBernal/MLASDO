% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performPCAAnalysis.R
\name{performPCAAnalysis}
\alias{performPCAAnalysis}
\title{performPCAAnalysis}
\usage{
performPCAAnalysis(
  savingName,
  mlAlgorithm,
  changedOmicData,
  selectedData,
  idColumn,
  classVariable
)
}
\arguments{
\item{savingName}{String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.}

\item{mlAlgorithm}{String | Machine Learning algorithm to be applied, the options are: Lasso or RF (Random Forest).}

\item{changedOmicData}{Dataset | Dataset of omic data that will be used.}

\item{selectedData}{Dataset | Dataset of omic data with only the predictors selected by the Lasso model.}

\item{idColumn}{String | Variable that indicates the identifier of each patient in both datasets. If the user does not specify a path to his own data, the value for the sample data, Trial, will be used.}

\item{classVariable}{String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.}
}
\description{
This function performs PCA analysis on the changed diagnoses after the execution of the genetic algorithm.
}
\examples{

MLASDO::performPCAAnalysis(savingName = savingName, mlAlgorithm = mlAlgorithm, changedOmicData = changedOmicData, selectedData = selectedData, idColumn = idColumn, classVariable = classVariable)
}
