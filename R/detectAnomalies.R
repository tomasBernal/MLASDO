#' executeAndPerformAnalysis
#'
#' Methodology based on the combination of a genetic algorithm and a Machine Learning technique to mutate the final
#' diagnosis of patients, detecting anomalies in them.
#'
#' These discrepancies may have their origin in the evolution of the subject itself, leading it from one group to
#' another in biological terms, or they may derive from human error in the labelling of the samples, among other causes.
#'
#' The identification and correction of these erroneous situations becomes essential to preserve the integrity
#' and accuracy of the data.
#'
#' @param justAnalysis Indicates whether to perform the analysis directly (TRUE) or to run the genetic algorithm (FALSE). Default value: FALSE.
#' @param solutionPath Path to genetic algorithm solution.
#'
#' @param mlAlgorithm Machine Learning algorithm to be applied, the options are: Lasso or RF (Random Forest). Default value: RF.
#'
#' @param numLassoExecutions Number of times the Lasso algorithm is executed. Default value: 5.
#' @param numTrees Number of trees of the Random Forest model. Default value: 100.
#' @param mtry Number of predictors that are evaluated at each partition (node) of each tree. Default value: 225.
#' @param splitrule This is the rule used by the algorithm to select the predictor and the optimal value to separate a node into two branches during the tree construction. Default value: gini.
#' @param sampleFraction Fraction of the training data that will be used to create each of the trees in the forest. Default value: 1.
#' @param maxDepth Maximum height of each tree in the forest. Default value: 4.
#' @param minNodeSize Minimum number of observations required in a node to be able to split it. Default value: 30.
#'
#' @param omicDataPath Path to omic data. If the user does not specify a path, the sample data will be used.
#' @param clinicDataPath Path to clinic data. If the user does not specify a path, the sample data will be used.
#' @param idColumn Variable that indicates the identifier of each patient in both datasets. If the user does not specify a path to his own data, the value for the sample data, Trial, will be used.
#' @param activePredictors Predictors on which the study of the ratios will be conducted after the genetic algorithm has been performed. Default value: All the predictors in clinic data, except classVariable and idColumn.
#' @param classVariable Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#' @param savingName Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#' @param nCores Number of cores to be used in parallelization. Default value: 6.
#' @param partitionPercentage Percentage (expressed as a fraction) with which the data will be split into a training and test set. Default value: 0.9 (90%).
#' @param nIterations Number of iterations (generations) the genetic algorithm will perform. Default value: 200.
#' @param nStopIter Number of iterations after which the algorithm will stop if all of them have the same fitness value. Default value: 25.
#' @param populationSize Number of solutions that will be part of the initial population. Default value: 150.
#' @param diagnosticChangeProbability Percentage (expressed as a fraction) indicating the probability of each gene in the solutions to be changed. Default value: 0.1 (10%).
#' @param crossoverOperator Crossover operator used in the genetic algorithm. Default value: Single Point Crossover.
#' @param crossoverProbability Percentage (expressed as a fraction) indicating the probability of crossover occurrence. Default value: 0.8 (80%).
#' @param selectionOperator Selection operator used in the genetic algorithm. Default value: Tournament Selection.
#' @param mutationOperator Mutation operator used in the genetic algorithm. Default value: Random Mutation.
#' @param mutationProbability Percentage (expressed as a fraction) indicating the probability of mutation occurrence. Default value: 0.1 (10%).
#' @param seed Seed used for the creation of training and test sets. Default value: 1234.
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::detectAnomalies(savingName = "DefaultExecutionLasso", mlAlgorithm = "Lasso",)
#'
#' MLASDO::detectAnomalies(savingName = "DefaultExecutionRF", mlAlgorithm = "RF")
#'
#' MLASDO::detectAnomalies(savingName = "QuickExecution", mlAlgorithm = "Lasso", nIterations = 3, nStopIter = 1, populationSize = 20, activePredictors = c("sex", "age", "Mutation", "Ethnicity"))
#'
#' MLASDO::detectAnomalies(savingName = "ExecutionWithOwnData", mlAlgorithm = "RF", omicDataPath = "./myOmicData.tsv", clinicDataPath = "./myClinicData.tsv", idColumn = "Patient.Id", nIterations = 3, populationSize = 10, classVariable = "Diagnosis", activePredictors = c("sex", "age", "Ethnicity"))
#'
#' MLASDO::detectAnomalies(justAnalysis = TRUE, solutionPath = "GA_solution.rds", savingName = "ExecutionWithOwnData", omicDataPath = "./myOmicData.tsv", clinicDataPath = "./myClinicData.tsv",idColumn = "Patient.Id",classVariable = "Diagnosis", activePredictors = c("sex", "age", "Ethnicity"))

detectAnomalies <- function(
    justAnalysis = FALSE,
    solutionPath = "",
    mlAlgorithm = "RF",
    numLassoExecutions = 5,
    numTrees = 100,
    mtry = 225,
    splitrule = "gini",
    sampleFraction = 1,
    maxDepth = 4,
    minNodeSize = 30,
    omicDataPath = "",
    clinicDataPath = "",
    idColumn,
    activePredictors = NULL,
    classVariable,
    savingName = "",
    nCores = 6,
    partitionPercentage = 0.9,
    nIterations = 200,
    nStopIter = 25,
    populationSize = 150,
    diagnosticChangeProbability = 0.1,
    crossoverOperator = gabin_spCrossover,
    crossoverProbability = 0.8,
    selectionOperator = gabin_tourSelection,
    mutationOperator = gabin_raMutation,
    mutationProbability = 0.1,
    seed = 1234
){


  #### CHECKING REQUIRED LIBRARIES ####
  if (!require("GA", character.only = TRUE)) {
    return("The package GA is not installed")
  }

  if (!require("ranger", character.only = TRUE)) {
    return("The package ranger is not installed")
  }

  if (!require("glmnet", character.only = TRUE)) {
    return("The package glmnet is not installed")
  }

  if (!require("caret", character.only = TRUE)) {
    return("The package caret is not installed")
  }

  if (!require("doParallel", character.only = TRUE)) {
    return("The package doParallel is not installed")
  }

  if (!require("dplyr", character.only = TRUE)) {
    return("The package dplyr is not installed")
  }

  #### CHECKING INITIAL STEP ####

  if(justAnalysis & solutionPath == ""){
    return("If you want to perform the analysis only, you must indicate the path to the final solution of the genetic algorithm.")
  }

  #### CHECKING ML ALGORITHM ####
  if(mlAlgorithm != "Lasso" & mlAlgorithm != "RF"){
    return("The Machine Learning algorithms that can be applied are: Lasso or RF (Random Forest)")
  }


  #### DATA READING ####

  if (omicDataPath == ""){
    omicRoute <- system.file("data", "omicData.tsv", package = "MLASDO")
    classVariable <- "Ca.Co.Last"
    idColumn <- "Trial"
  } else {
    omicRoute <- omicDataPath
  }

  if (clinicDataPath == ""){
    clinicRoute <- system.file("data", "clinicData.tsv", package = "MLASDO")
    classVariable <- "Ca.Co.Last"
    idColumn <- "Trial"
  } else {
    clinicRoute <- clinicDataPath
  }

  omicData <- read.table(omicRoute, header = TRUE, sep = "\t", row.names = 1)
  clinicData <- read.table(clinicRoute, header = TRUE, sep = "\t", row.names = 1)

  #### CHECKING PARAMETERS ####

  if(is.null(activePredictors)){
    activePredictors <- names(clinicData)

    activePredictors <- activePredictors[activePredictors != classVariable]

    activePredictors <- activePredictors[activePredictors != idColumn]
  }

  # Checking if the patients in both datasets are the same
  omicIds <- sort(omicData[[idColumn]])
  clinicIDs <- sort(clinicData[[idColumn]])

  if(!identical(omicIds, clinicIDs)){
    return("The IDs of the patients don't match in both datasets")
  }

  # First, check if the selected variable as a predictor is included in the variables to be studied
  if (classVariable %in% activePredictors) {
    return("The selected variable for prediction cannot be included in the variables chosen for analysis")
  }

  if (idColumn %in% activePredictors) {
    return("The variable that identifies the patients cannot be included in the variables chosen for analysis")
  }

  # After this, check if the selected variable as the class variable for prediction is valid
  if (!classVariable %in% names(omicData)) {
    return("The selected variable for prediction does not exist in the omic dataframe")
  }

  if (!classVariable %in% names(clinicData)) {
    return("The selected variable for prediction does not exist in the clinic dataframe")
  }

  if (!idColumn %in% names(clinicData)) {
    return("The variable that identifies the patients does not exist in the clinic dataframe")
  }

  if (!idColumn %in% names(omicData)) {
    return("The variable that identifies the patients does not exist in the clinic omic dataframe")
  }

  # Next, check if the selected variable as the class variable for prediction is binary
  if (length(unique(omicData[[classVariable]])) > 2) {
    return("The selected variable for prediction is not binary")
  }

  if (length(unique(clinicData[[classVariable]])) > 2) {
    return("The selected variable for prediction is not binary")
  }

  # Checking if the class variable is well coded
  if(all(!(omicData[[classVariable]] %in% c(0, 1))) & all(!(omicData[[classVariable]] %in% c("Case", "Control")))){
    return("The class variable is not well coded")
  }

  if(all(!(clinicData[[classVariable]] %in% c(0, 1))) & all(!(clinicData[[classVariable]] %in% c("Case", "Control")))){
    return("The class variable is not well coded")
  }

  # Checking if the coding of the class variable in the omic dataset is correct
  if (0 %in% unique(omicData[[classVariable]])) {

    print("The class variable is coded with 0, those values are changed to Control.")

    omicData[[classVariable]][omicData[[classVariable]] == 0] <- "Control"

    clinicData[[classVariable]][clinicData[[classVariable]] == 0] <- "Control"

  }

  # Checking if the coding of the class variable in the omic dataset is correct
  if (1 %in% unique(omicData[[classVariable]])) {

    print("The class variable is coded with 1, those values are changed to Case.")

    omicData[[classVariable]][omicData[[classVariable]] == 1] <- "Case"

    clinicData[[classVariable]][clinicData[[classVariable]] == 1] <- "Case"

  }

  # Checking if the coding of the class variable in the clinic dataset is correct
  if (0 %in% unique(clinicData[[classVariable]])) {

    print("The class variable is coded with 0, those values are changed to Control.")

    omicData[[classVariable]][omicData[[classVariable]] == 0] <- "Control"

    clinicData[[classVariable]][clinicData[[classVariable]] == 0] <- "Control"

  }

  # Checking if the coding of the class variable in the clinic dataset is correct
  if (1 %in% unique(clinicData[[classVariable]])) {

    print("The class variable is coded with 1, those values are changed to Case.")

    omicData[[classVariable]][omicData[[classVariable]] == 1] <- "Case"

    clinicData[[classVariable]][clinicData[[classVariable]] == 1] <- "Case"

  }

  # Finally, check if the variables selected for ratio analysis are valid
  invalidPredictors <- activePredictors[!(activePredictors %in% names(clinicData))]

  if (length(invalidPredictors) > 0) {
    cat("The following selected variables for analysis do not exist in the dataframe:\n")
    return(paste(invalidPredictors, collapse = ", "))
  }

  # If the name is empty, we need to establish a default one
  if (savingName == ""){
    savingName <- format(Sys.time(), "%d-%m-%H-%M")
  }


  #### DATA PROCESSING ####

  omicData <- omicData[,!(names(omicData) %in% activePredictors)]

  omicsFiltered <- omicData[,!(names(omicData) %in% classVariable)]

  validClinicData <- merge(clinicData, omicsFiltered, by = idColumn)

  validClinicData <- validClinicData[, (names(validClinicData) %in% c(activePredictors, idColumn, classVariable))]

  if(!justAnalysis){
    print("Executing the genetic algorithm")
    MLASDO::executeGA(mlAlgorithm = mlAlgorithm, numLassoExecutions = numLassoExecutions, numTrees = numTrees, mtry = mtry, splitrule = splitrule, sampleFraction = sampleFraction, maxDepth = maxDepth, minNodeSize = minNodeSize, omicData = omicData, savingName = savingName, classVariable = classVariable, activePredictors = activePredictors, nCores = nCores, partitionPercentage = partitionPercentage, nIterations = nIterations, nStopIter = nStopIter, populationSize = populationSize, diagnosticChangeProbability = diagnosticChangeProbability, crossoverOperator = crossoverOperator, crossoverProbability = crossoverProbability, selectionOperator = selectionOperator, mutationOperator = mutationOperator, mutationProbability = mutationProbability, seed = seed)
  }

  print("Performing PCA analysis")
  MLASDO::performPCAAnalysis(justAnalysis = justAnalysis, solutionPath = solutionPath, idColumn = idColumn, omicData = omicData, savingName = savingName, classVariable = classVariable, activePredictors = activePredictors)

  print("Performing ratio analysis")
  MLASDO::performRatioAnalysis(justAnalysis = justAnalysis, solutionPath = solutionPath, clinicData = validClinicData, savingName = savingName, classVariable = classVariable, activePredictors = activePredictors)


  # Reading the GA solution
  name <- paste("GA", savingName, sep="_")

  if(justAnalysis){
    solutionGA <- readRDS(solutionPath)
  } else {
    gaPath <- paste(name, "Solution.rds", sep="_")
    solutionGA <- readRDS(gaPath)
  }


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

  changedClinicData <- validClinicData

  changedClinicData[[classVariable]] <- changedDiagnoses

  print("Compiling Markdown file")
  MLASDO::compileMarkdown(savingName = savingName, clinicData = changedClinicData, classVariable = classVariable)
}
