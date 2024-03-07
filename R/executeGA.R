#' executeGA
#'
#' @description Methodology based on the combination of a genetic algorithm and a Machine Learning technique to mutate the final
#' diagnosis of patients, detecting anomalies in them.
#'
#' These discrepancies may have their origin in the evolution of the subject itself, leading it from one group to
#' another in biological terms, or they may derive from human error in the labelling of the samples, among other causes.
#'
#' This detection is done by combining a genetic algorithm and a machine learning algorithm, it is possible to use a Lasso model
#' or a Random Forest model.
#'
#' A genetic algorithm is a technique inspired by biological evolution that simulates the process of natural selection and
#' evolution of species, the execution of a genetic algorithm starts with a randomly generated initial solution set.

#' - In this case, it is assumed that certain patients are misdiagnosed, so the initial solution set will be created by
#' randomly changing 10% of the patient diagnoses in each of these solutions.

#' In each generation (iteration) of the genetic algorithm, each of the solutions from the set of possible solutions is
#' used as a training dataset for the ML model, which will classify the patients and whose performance will determine the
#' fitness value (score) of the solution.

#' The way in which the score for each solution is calculated depends on the model used:

#'  - For the Lasso model, the score of each solution is defined by the average of the balanced accuracies obtained after
#'  training 5 Lasso models with that solution as training dataset.

#' - For RF model, the score of each solution is defined by the balanced accuracy obtained after training the RF model with
#' that solution as training dataset.

#' After this evaluation, a ranking of the solutions is created:

#'   - The two best solutions are combined to generate a new solution that is added to the set of possible solutions to the
#'   problem.

#' - Simultaneously, a certain number of the worst solutions are eliminated from the generation, i.e. they are no longer part
#' of the set of possible solutions to the problem.

#' This process is repeated for the specified number of generations, until a certain predefined threshold is reached in the
#' score of the solutions or until the value of the best score of the set of solutions is repeated for a certain specified
#' number of generations.

#' This procedure is used to determine which samples were anomalous. The genetic algorithm locates these anomalous samples
#' and the machine learning algorithm evidences this detection, validating that these samples belong to the alternative group.
#'
#'
#' @param mlAlgorithm String | Machine Learning algorithm to be applied, the options are: Lasso or RF (Random Forest). Default value: RF.
#'
#' @param numModelExecutions Integer | Number of times the Lasso algorithm is executed. Default value: 5.
#'
#' @param predictorsToSelect Integer | Number of predictors to be selected from the most important predictors ranked by the RF model. This parameter is a integer number between 1 and the total number of predictors in the data. Default value: 15.
#' @param numTrees Integer | Number of trees of the Random Forest model. Default value: 100.
#' @param mtry Integer | Number of predictors that are evaluated at each partition (node) of each tree. Default value: 225.
#' @param splitrule String | This is the rule used by the algorithm to select the predictor and the optimal value to separate a node into two branches during the tree construction. Default value: gini.
#' @param sampleFraction Decimal | Fraction of the training data that will be used to create each of the trees in the forest. Default value: 1.
#' @param maxDepth Integer | Maximum height of each tree in the forest. Default value: 4.
#' @param minNodeSize Integer | Minimum number of observations required in a node to be able to split it. Default value: 30.
#'
#' @param omicData Dataset of omic data that will be used.
#' @param subsetTrain Subset of train samples.
#' @param idColumn String | Variable that indicates the identifier of each patient in both datasets. If the user does not specify a path to his own data, the value for the sample data, Trial, will be used.
#' @param activePredictors Array of Strings | Predictors on which the study of the ratios will be conducted after the genetic algorithm has been performed. Default value: All the predictors in clinic data, except classVariable and idColumn.
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#' @param nCores Integer | Number of cores to be used in parallelization. Default value: 6.
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
#' MLASDO::executeGA(mlAlgorithm = mlAlgorithm, numModelExecutions = numModelExecutions, predictorsToSelect = predictorsToSelect, numTrees = numTrees, mtry = mtry, splitrule = splitrule, sampleFraction = sampleFraction, maxDepth = maxDepth, minNodeSize = minNodeSize, omicData = omicData, subsetTrain = subsetTrain, activePredictors = activePredictors, classVariable = classVariable, idColumn = idColumn, savingName = savingName, nCores = nCores, nIterations = nIterations, nStopIter = nStopIter, populationSize = populationSize, diagnosticChangeProbability = diagnosticChangeProbability, crossoverOperator = crossoverOperator, crossoverProbability = crossoverProbability, selectionOperator = selectionOperator, mutationOperator = mutationOperator, mutationProbability = mutationProbability, seed = seed)

executeGA <- function(
                      mlAlgorithm,
                      numModelExecutions,
                      predictorsToSelect,
                      numTrees,
                      mtry,
                      splitrule,
                      sampleFraction,
                      maxDepth,
                      minNodeSize,
                      omicData,
                      subsetTrain,
                      activePredictors,
                      classVariable,
                      idColumn,
                      savingName,
                      nCores,
                      nIterations,
                      nStopIter,
                      populationSize,
                      diagnosticChangeProbability,
                      crossoverOperator,
                      crossoverProbability,
                      selectionOperator,
                      mutationOperator,
                      mutationProbability,
                      seed
                      ){

  #### REQUIRED LIBRARIES ####
  library(GA) # For genetic algorithm
  library(glmnet) # For Lasso model
  library(ranger) # For ranger model
  library(caret) # For confusion matrix
  library(doParallel) # For parallel execution


  #### DATA READING ####
  omic <- omicData

  # We need to remove the column that indicates the patient id
  omic[[idColumn]] <- NULL

  omic[[classVariable]] <- ifelse(omic[[classVariable]] == "Case", 1, 0)

  #### CREATION OF TRAIN AND TEST SETS ####
  omicTrain <- omic[subsetTrain,]

  omicTrainDiagnosis <- omicTrain[[classVariable]]
  omicTrain[[classVariable]] <- NULL

  omicTest <- omic[-subsetTrain,]

  omicTestDiagnosis <- omicTest[[classVariable]]
  omicTest[[classVariable]] <- NULL


  #### GENETIC ALGORITHM CREATION AND EXECUTION ####

  ## Definition of the fitness function
  # This function evaluates each possible solution
  fitness <- function(genome) {

    # Substitute the NA values for 0
    genome <- replace(genome, is.na(genome), 0)

    # Calculating new diagnoses determined by the current solution
    # It is done through an XOR operation
    solutionData <- bitwXor(omicTrainDiagnosis, genome)

    # This vector will store the balanced means of the numModelExecutions executions
    balancedAccValues <- vector(length=numModelExecutions)

    # Train the Lasso model numModelExecutions times
    for (i in 1:numModelExecutions) {

        if(mlAlgorithm == "Lasso"){

          set.seed(seed)

          # Train the Lasso model
          model <- cv.glmnet(as.matrix(omicTrain), as.matrix(solutionData), alpha = 1, family = "binomial", type.measure = "class", nfolds = 10)

          # Use the model to predict on the test set
          modelPrediction <- predict(model, newx = as.matrix(omicTest), alpha = 1, s = "lambda.min", type = "class")

        } else if(mlAlgorithm == "RF"){

          # Creating the ranger model
          model <- ranger(
            x = omicTrain,
            y = solutionData,
            num.trees = numTrees,
            mtry = mtry,
            splitrule = splitrule,
            importance = "impurity", # In order to obtain the important variables in the prediction
            sample.fraction = sampleFraction,
            max.depth = maxDepth,
            min.node.size = minNodeSize,
            classification = TRUE,
            seed = seed
          )

          # Use the model to predict on the test set
          modelPrediction <- predict(model, omicTest)$predictions
        }

        # Get the confusion matrix
        cfModel <- confusionMatrix(as.factor(as.integer(modelPrediction)), as.factor(omicTestDiagnosis))

        specificity <- ifelse(is.na(as.numeric(cfModel$byClass["Specificity"])), 0,  as.numeric(cfModel$byClass["Specificity"]))

        sensitivity <- ifelse(is.na(as.numeric(cfModel$byClass["Sensitivity"])), 0,  as.numeric(cfModel$byClass["Sensitivity"]))

        # Save the balanced mean obtained in this iteration
        balancedAccValues[i] <- (specificity + sensitivity) / 2

      }

      # Return the balanced accuracy of the ranger model
      return(mean(balancedAccValues))

  }


  ## Definition of the function that generates initial population
  # Function to generate the initial population of individuals
  generateInitialPopulation <- function(object) {

    population <- matrix(
                      as.double(NA),
                      nrow = object@popSize,
                      ncol = object@nBits
                  )

    for(j in 1:object@popSize){

      population[j, ] <- round(rbinom(object@nBits, size = 1, prob = diagnosticChangeProbability))
    }

    storage.mode(population) <- "integer"
    return(population)
  }


  ## Genetic Algorithm Definition
  GA <- ga(
    keepBest = TRUE, # To save the best solutions from each iteration
    type = "binary", # With genome encoding based on binary (0 = Maintain diagnosis, 1 = Change diagnosis)
    maxFitness = 1, # Maximum fitness value (to stop if reached)
    fitness = fitness, # Fitness function (explained above)
    seed = seed, # Seed used by the genetic algorithm in random processes (crossover and mutation)
    nBits = nrow(omicTrain), # Number of bits in the genome (one for each example in the dataset)
    popSize = populationSize, # Size of the initial population
    population = generateInitialPopulation, # Function that generates the initial population (explained above)
    selection = selectionOperator, # Selection of the best individuals through tournament
    maxiter = nIterations, # Number of iterations to perform
    run = nStopIter, # Number of iterations to stop if the fitness value doesn't change
    pcrossover = crossoverProbability, # Probability of performing crossover between parents
    crossover = crossoverOperator, # Function that performs crossover, in this case, single-point crossover (taking a part from each parent)
    pmutation = mutationProbability, # Probability of performing mutation in parents
    mutation = mutationOperator, # Function that performs mutation, in this case, randomly changes a gene
    parallel = nCores # Number of cores used in the parallelization
  )


  ## Save the result of the genetic algorithm
  dirPath <- paste(savingName, "geneticAlgorithm", sep = "/")

  name <- paste("GA", savingName, sep="_")
  gaPath <- paste(name, ".rds", sep="")
  saveRDS(GA, file = paste(dirPath, gaPath, sep = "/"))

  ## Obtain the best solution obtained by the genetic algorithm and its accuracy
  maxValue <- 0

  # Parallelizing
  cl <- makeCluster(nCores)
  registerDoParallel(cl)

  for(i in 1:nrow(GA@solution)){

    currentVal <- fitness((GA@solution[i,]))

    if(currentVal > maxValue){
      geneticSolution <- GA@solution[i,]
      maxValue <- currentVal
    }
  }

  # Stop parallelization
  stopCluster(cl)

  ## Save the best solution of the genetic algorithm
  solutionPath <- paste(name, "Solution.rds", sep="_")
  saveRDS(geneticSolution, file = paste(dirPath, solutionPath, sep = "/"))

  # Obtaining the worst and the best model with the final solution
  solutionData <- bitwXor(omicTrainDiagnosis, geneticSolution)

  models <- vector(mode = "list", length = numModelExecutions)

  predictorsInfo <- list()

  bestModelIndex <- 1
  worstModelIndex <- 1

  bestModelBA <- 0
  worstModelBA <- 1


  for (i in 1:numModelExecutions) {

    if(mlAlgorithm == "Lasso"){

      set.seed(seed)

      # Train the Lasso model
      model <- cv.glmnet(as.matrix(omicTrain), as.matrix(solutionData), alpha = 1, family = "binomial", type.measure = "class", nfolds = 10)

      # Use the model to predict on the test set
      modelPrediction <- predict(model, newx = as.matrix(omicTest), alpha = 1, s = "lambda.min", type = "class")

      predInfo <- coef(model, s = model$lambda.min)

      indexes <- predInfo@i

      posZero <- which(indexes == 0)

      indexes <- indexes[-posZero]

      for(i in 1:length(indexes)){

        predName <- names(omicTrain)[[indexes[[i]]]]

        predImp <- predInfo@x[[i]]

        if (exists(predName, where = predictorsInfo)) {

          actualAparitions <- predictorsInfo[[predName]][1]

          importance <- predictorsInfo[[predName]][2:length(predictorsInfo[[predName]])]

          predictorsInfo[[predName]] <- c(actualAparitions + 1, importance, predImp)

        } else {

          predictorsInfo[[predName]] <- c(1, predImp)

        }

      }

    } else if(mlAlgorithm == "RF"){

      # Creating the ranger model
      model <- ranger(
        x = omicTrain,
        y = solutionData,
        num.trees = numTrees,
        mtry = mtry,
        splitrule = splitrule,
        importance = "impurity", # In order to obtain the important variables in the prediction
        sample.fraction = sampleFraction,
        max.depth = maxDepth,
        min.node.size = minNodeSize,
        classification = TRUE,
        seed = seed
      )

      # Use the model to predict on the test set
      modelPrediction <- predict(model, omicTest)$predictions

      varImportance <- model$variable.importance

      predInfo <- data.frame(Variable = names(varImportance), Importance = as.numeric(varImportance))

      predInfo <- predInfo[order(-predInfo$Importance), ]

      predInfo <- head(predInfo, predictorsToSelect)

      predInfo <- predInfo[order(predInfo$Importance), ]

      for(i in 1:predictorsToSelect){

        if (exists(predInfo$Variable[[i]], where = predictorsInfo)) {

          actualAparitions <- predictorsInfo[[predInfo$Variable[[i]]]][1]

          importance <- predictorsInfo[[predInfo$Variable[[i]]]][2:length(predictorsInfo[[predInfo$Variable[[i]]]])]

          predictorsInfo[[predInfo$Variable[[i]]]] <- c(actualAparitions + 1, importance, abs(predInfo$Importance[[i]]))

        } else {

          predictorsInfo[[predInfo$Variable[[i]]]] <- c(1, abs(predInfo$Importance[[i]]))

        }
      }
    }

    models[[i]] <- model

    # Get the confusion matrix
    cfModel <- confusionMatrix(as.factor(as.integer(modelPrediction)), as.factor(omicTestDiagnosis))

    specificity <- ifelse(is.na(as.numeric(cfModel$byClass["Specificity"])), 0,  as.numeric(cfModel$byClass["Specificity"]))

    sensitivity <- ifelse(is.na(as.numeric(cfModel$byClass["Sensitivity"])), 0,  as.numeric(cfModel$byClass["Sensitivity"]))

    # Save the balanced mean obtained in this iteration
    actualBA <- (specificity + sensitivity) / 2

    if(actualBA > bestModelBA){

      bestModelBA <- actualBA
      bestModelIndex <- i
    }

    if(actualBA < worstModelBA){
      worstModelBA <- actualBA
      worstModelIndex <- i
    }

  }


  ## Save the best solution of the genetic algorithm
  modelPath <- paste(name, "Best_Model.rds", sep="_")
  saveRDS(models[[bestModelIndex]], file = paste(dirPath, modelPath, sep = "/"))

  modelPath <- paste(name, "Worst_Model.rds", sep="_")
  saveRDS(models[[worstModelIndex]], file = paste(dirPath, modelPath, sep = "/"))


  predictorsInfo <- as.data.frame(predictorsInfo)

  modelCols <- vector("character", numModelExecutions)

  for(i in 1:numModelExecutions){

    modelCols[[i]] <- paste("model", i, sep = "_")
  }

  newColnames <- c("numAparitions", modelCols)

  predictorsInfo <- t(predictorsInfo)

  colnames(predictorsInfo) <- newColnames

  meanImportance <- vector(mode = "list", length = nrow(predictorsInfo))

  for(i in 1:nrow(predictorsInfo)){

    meanImportance[[i]] <- mean(predictorsInfo[i, 2:ncol(predictorsInfo)])
  }

  predictorsImportancePath <- paste(name, "Predictors_Importance_a.tsv", sep="_")

  write.table(predictorsInfo, paste(dirPath, predictorsImportancePath, sep = "/"), row.names = T, col.names = T, sep =  '\t')

  print(unlist(meanImportance))

  predictorsInfo$meanImportance <- unlist(meanImportance)

  predictorsImportancePath <- paste(name, "Predictors_Importance.tsv", sep="_")

  write.table(predictorsInfo, paste(dirPath, predictorsImportancePath, sep = "/"), row.names = T, col.names = T, sep =  '\t')

  postFitness <- function(genome) {

    # Substitute the NA values for 0
    genome <- replace(genome, is.na(genome), 0)

    # Calculating new diagnoses determined by the current solution
    # It is done through an XOR operation
    solutionData <- bitwXor(omicTrainDiagnosis, genome)

    # This vector will store the balanced means of the numModelExecutions executions
    balancedAccValues <- vector(length = numModelExecutions)
    numPredictors <- vector(length = numModelExecutions)

    set.seed(seed)

    # Train the Lasso model numModelExecutions times
    for (i in 1:numModelExecutions) {

      # Train the Lasso model
      model <- cv.glmnet(as.matrix(omicTrain), as.matrix(solutionData), alpha = 1, family = "binomial", type.measure = "class", nfolds = 10)

      coefficients <- coef(model, s = model$lambda.min)

      numPredictors[i] <- length(coefficients@i)

      # Use the model to predict on the test set
      modelPrediction <- predict(model, newx = as.matrix(omicTest), alpha = 1, s = "lambda.min", type = "class")

      # Get the confusion matrix
      cfModel <- confusionMatrix(as.factor(as.integer(modelPrediction)), as.factor(omicTestDiagnosis))

      specificity <- ifelse(is.na(as.numeric(cfModel$byClass["Specificity"])), 0,  as.numeric(cfModel$byClass["Specificity"]))

      sensitivity <- ifelse(is.na(as.numeric(cfModel$byClass["Sensitivity"])), 0,  as.numeric(cfModel$byClass["Sensitivity"]))

      # Save the balanced mean obtained in this iteration
      balancedAccValues[i] <- (specificity + sensitivity) / 2

    }

    # Return the mean of the numModelExecutions balanced means obtained
    return(c(mean(balancedAccValues), round(mean(numPredictors))))
  }

  if(mlAlgorithm == "Lasso"){

    # Parallelizing
    cl <- makeCluster(nCores)
    registerDoParallel(cl)

    numberPredictorsLasso <- c()

    for(i in 1:length(GA@bestSol)){

      maxValue <- 0
      maxPredictors <- 0

      for(j in 1:nrow(GA@bestSol[[i]])){

        currentVal <- postFitness((GA@bestSol[[i]][j ,]))

        if(currentVal[[1]] > maxValue){
          maxValue <- currentVal[[1]]
          maxPredictors <- currentVal[[2]]
        }

      }

      numberPredictorsLasso <- c(numberPredictorsLasso, maxPredictors)

    }

    # Stop parallelization
    stopCluster(cl)

    ## Save the best solution of the genetic algorithm
    predictorsPath <- paste(name, "Predictors.rds", sep="_")
    saveRDS(numberPredictorsLasso, file = paste(dirPath, predictorsPath, sep = "/"))

  }
}
