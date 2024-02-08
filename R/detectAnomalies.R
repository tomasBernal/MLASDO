#' executeAndPerformAnalysis
#'
#' @description Methodology based on the combination of a genetic algorithm and a Machine Learning technique to mutate the final
#' diagnosis of patients, detecting anomalies in them.
#'
#' These discrepancies may have their origin in the evolution of the subject itself, leading it from one group to
#' another in biological terms, or they may derive from human error in the labelling of the samples, among other causes.
#'
#' The identification and correction of these erroneous situations becomes essential to preserve the integrity
#' and accuracy of the data.
#'
#' @param justAnalysis Bool | Indicates whether to perform the analysis directly (TRUE) or to run the genetic algorithm (FALSE). Default value: FALSE.
#' @param geneticPath String | Path to genetic algorithm object.
#' @param solutionPath String | Path to genetic algorithm solution.
#' @param modelPath String | Path to the best model obtained.
#' @param lassoPredictorsPath String | Path to the mean number of predictors selected by Lasso in each generation.
#'
#' @param mlAlgorithm String | Machine Learning algorithm to be applied, the options are: Lasso or RF (Random Forest).
#'
#' @param numLassoExecutions Integer | Number of times the Lasso algorithm is executed. Default value: 5.
#' @param numTrees Integer | Number of trees of the Random Forest model. Default value: 100.
#' @param mtry Integer | Number of predictors that are evaluated at each partition (node) of each tree. Default value: 225.
#' @param splitrule String | This is the rule used by the algorithm to select the predictor and the optimal value to separate a node into two branches during the tree construction. Default value: gini.
#' @param sampleFraction Decimal | Fraction of the training data that will be used to create each of the trees in the forest. Default value: 1.
#' @param maxDepth Integer | Maximum height of each tree in the forest. Default value: 4.
#' @param minNodeSize Integer | Minimum number of observations required in a node to be able to split it. Default value: 30.
#'
#' @param omicDataPath String | Path to omic data. If the user does not specify a path, the sample data will be used.
#' @param clinicDataPath String | Path to clinic data. If the user does not specify a path, the sample data will be used.
#' @param omicPredictorsToIgnore Array of Strings | Variables to be removed from the omic dataset. These will not be taken into account in the execution.
#' @param clinicPredictorsToIgnore Array of Strings | Variables to be removed from the clinic dataset. These will not be taken into account in the execution.
#'
#' @param idColumn String | Variable that indicates the identifier of each patient in both datasets. If the user does not specify a path to his own data, the value for the sample data, Trial, will be used.
#' @param activePredictors Array of Strings | Predictors on which the study of the ratios will be conducted after the genetic algorithm has been performed. Default value: All the predictors in clinic data, except classVariable and idColumn.
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#' @param nCores Integer | Number of cores to be used in parallelization. Default value: 6.
#' @param partitionPercentage Decimal | Percentage (expressed as a fraction) with which the data will be split into a training and test set. Default value: 0.9 (90%).
#' @param nIterations Integer | Number of iterations (generations) the genetic algorithm will perform. Default value: 200.
#' @param nStopIter Integer | Number of iterations after which the algorithm will stop if all of them have the same fitness value. Default value: 25.
#' @param populationSize Integer | Number of solutions that will be part of the initial population. Default value: 150.
#' @param diagnosticChangeProbability Decimal | Percentage (expressed as a fraction) indicating the probability of each gene in the solutions to be changed. Default value: 0.1 (10%).
#' @param crossoverOperator String | Crossover operator used in the genetic algorithm. Default value: Single Point Crossover.
#' @param crossoverProbability Decimal | Percentage (expressed as a fraction) indicating the probability of crossover occurrence. Default value: 0.8 (80%).
#' @param selectionOperator String | Selection operator used in the genetic algorithm. Default value: Tournament Selection.
#' @param mutationOperator String | Mutation operator used in the genetic algorithm. Default value: Random Mutation.
#' @param mutationProbability Decimal | Percentage (expressed as a fraction) indicating the probability of mutation occurrence. Default value: 0.1 (10%).
#' @param seed Integer | Seed used for the creation of training and test sets. Default value: 1234.
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
#' MLASDO::detectAnomalies(savingName = "QuickExecution", mlAlgorithm = "Lasso", nIterations = 3, nStopIter = 2, populationSize = 20, activePredictors = c("sex", "age", "Mutation", "Ethnicity"))
#'
#' MLASDO::detectAnomalies(savingName = "ExecutionWithOwnData", mlAlgorithm = "RF", omicDataPath = "./myOmicData.tsv", clinicDataPath = "./myClinicData.tsv", idColumn = "Patient.Id", nIterations = 3, populationSize = 10, classVariable = "Diagnosis", activePredictors = c("sex", "age", "Ethnicity"))
#'
#' MLASDO::detectAnomalies(justAnalysis = TRUE, mlAlgorithm = "Lasso", geneticPath = "GA.rds", solutionPath = "GA_solution.rds", modelPath = "GA_Model.rds", lassoPredictorsPath = "GA_Lasso_Predictors.rds", savingName = "ExecutionWithOwnData", omicDataPath = "./myOmicData.tsv", clinicDataPath = "./myClinicData.tsv",idColumn = "Patient.Id",classVariable = "Diagnosis", activePredictors = c("sex", "age", "Ethnicity"))

detectAnomalies <- function(
    justAnalysis = FALSE,
    solutionPath = "",
    geneticPath = "",
    modelPath = "",
    lassoPredictorsPath = "",
    mlAlgorithm,
    numLassoExecutions = 5,
    numTrees = 100,
    mtry = 225,
    splitrule = "gini",
    sampleFraction = 1,
    maxDepth = 4,
    minNodeSize = 30,
    omicDataPath = "",
    clinicDataPath = "",
    omicPredictorsToIgnore = NULL,
    clinicPredictorsToIgnore = NULL,
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

  if (!require("plotly", character.only = TRUE)) {
    return("The package plotly is not installed")
  }

  if (!require("lazyeval", character.only = TRUE)) {
    return("The package lazyeval is not installed")
  }

  #### CHECKING INITIAL STEP ####

  if(justAnalysis & solutionPath == ""){
    return("If you want to perform the analysis only, you must indicate the path to the final solution of the genetic algorithm.")
  }

  if(justAnalysis & geneticPath == ""){
    return("If you want to perform the analysis only, you must indicate the path to the genetic algorithm object.")
  }

  if(justAnalysis & modelPath == ""){
    return("If you want to perform the analysis only, you must indicate the path to the model.")
  }

  if(justAnalysis & mlAlgorithm == "Lasso" & lassoPredictorsPath == ""){
    return("If you want to perform the analysis of a Lasso model, you must indicate the path to the mean number of predictors selected in each generation.")
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

  if(!is.null(omicPredictorsToIgnore)){

    validColumnsToRemove <- names(omicData) %in% omicPredictorsToIgnore

    if(any(validColumnsToRemove)){
      omicData <- omicData[, !validColumnsToRemove]
    }

  }

  if(!is.null(clinicPredictorsToIgnore)){

    validColumnsToRemove <- names(clinicData) %in% clinicPredictorsToIgnore

    if(any(validColumnsToRemove)){
      clinicData <- clinicData[, !validColumnsToRemove]
    }

  }

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
    print("The following selected variables for analysis do not exist in the dataframe:\n")
    print(paste(invalidPredictors, collapse = ", "))

    activePredictors <- activePredictors[!(activePredictors %in% invalidPredictors)]
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

  dir.create(savingName)
  dir.create(paste(savingName, "analysisData", sep = "/"))

  set.seed(seed)
  subsetTrain <- sample(1:nrow(omicData), nrow(omicData) * partitionPercentage)

  if(!justAnalysis){

    dir.create(paste(savingName, "geneticAlgorithm", sep = "/"))

    print("Executing the genetic algorithm")
    MLASDO::executeGA(
      mlAlgorithm = mlAlgorithm,
      numLassoExecutions = numLassoExecutions,
      numTrees = numTrees,
      mtry = mtry,
      splitrule = splitrule,
      sampleFraction = sampleFraction,
      maxDepth = maxDepth,
      minNodeSize = minNodeSize,
      omicData = omicData,
      subsetTrain = subsetTrain,
      savingName = savingName,
      classVariable = classVariable,
      activePredictors = activePredictors,
      nCores = nCores,
      nIterations = nIterations,
      nStopIter = nStopIter,
      populationSize = populationSize,
      diagnosticChangeProbability = diagnosticChangeProbability,
      crossoverOperator = crossoverOperator,
      crossoverProbability = crossoverProbability,
      selectionOperator = selectionOperator,
      mutationOperator = mutationOperator,
      mutationProbability = mutationProbability,
      seed = seed
      )
  }

  # Reading the GA solution
  name <- paste("GA", savingName, sep="_")

  if(justAnalysis){

    geneticAlgorithm <- readRDS(geneticPath)
    solutionGA <- readRDS(solutionPath)
    model <- readRDS(modelPath)

  } else {

    dirPath <- paste(savingName, "geneticAlgorithm", name, sep = "/")

    modelPath <- paste(dirPath, "Model.rds", sep="_")
    model <- readRDS(modelPath)

    solutionPath <- paste(dirPath, "Solution.rds", sep="_")
    solutionGA <- readRDS(solutionPath)

    gaPath <- paste(dirPath, ".rds", sep="")
    geneticAlgorithm <- readRDS(gaPath)

  }

  # Creating a copy of the original diagnosis
  changedDiagnoses <- omicData[subsetTrain,][[classVariable]]

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

  changedOmicData <- omicData

  changedClinicData <- validClinicData

  changedClinicData[subsetTrain,][[classVariable]] <- changedDiagnoses
  changedOmicData[subsetTrain,][[classVariable]] <- changedDiagnoses

  print("Performing PCA analysis")
  MLASDO::performPCAAnalysis(
    model = model,
    idColumn = idColumn,
    mlAlgorithm = mlAlgorithm,
    changedOmicData = changedOmicData,
    savingName = savingName,
    classVariable = classVariable
  )

  print("Performing ratio analysis")
  MLASDO::performRatioAnalysis(
    changedClinicData = changedClinicData,
    firstGroup = firstGroup,
    secondGroup = secondGroup,
    savingName = savingName,
    classVariable = classVariable,
    activePredictors = activePredictors
  )

  print("Compiling Markdown file")
  MLASDO::compileMarkdown(
    savingName = savingName,
    justAnalysis = justAnalysis,
    lassoPredictorsPath = lassoPredictorsPath,
    model = model,
    mlAlgorithm = mlAlgorithm,
    geneticAlgorithm = geneticAlgorithm,
    originalDiagnosis = omicData[[classVariable]],
    clinicData = changedClinicData,
    classVariable = classVariable
    )
}
