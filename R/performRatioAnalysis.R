#' performRatioAnalysis
#'
#' This function performs ratio-based analysis on the active activePredictors passed as parameters.
#'
#'
#' @param clinicData Dataset of clinic data that will be used.
#' @param activePredictors activePredictors on which the study of the ratios will be conducted after the genetic algorithm has been performed.
#' @param classVariable Target variable, which must be binary, meaning it has two possible values.
#' @param savingName Name under which the model and solution were saved after the execution of the genetic algorithm.
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::performRatioAnalysis(clinicData = clinicData, activePredictors = activePredictors, classVariable = classVariable, savingName = savingName)


performRatioAnalysis <- function(
    clinicData,
    activePredictors,
    classVariable,
    savingName
){

  #### DATA READING ####
  clinic <- clinicData

  #### GENETIC ALGORITHM SOLUTION READING ####
  name <- paste("GA", savingName, sep="_")
  solutionPath <- paste(name, "Solution.rds", sep="_")
  solutionGA <- readRDS(solutionPath)

  #### DATA PROCESSING ####

  # Creating a copy of the original diagnosis
  changedDiagnoses <- clinic[[classVariable]]

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

  # Saving the changes
  clinic[[classVariable]] <- changedDiagnoses

  # Iterating through the selected activePredictors for analysis
  for (activePredictor in activePredictors){

    # Checking if the activePredictor is categorical
    if(!is.numeric(clinic[[activePredictor]])){


        # Replacing NA with the string "Empty"
        clinic[[activePredictor]] <- replace(clinic[[activePredictor]], is.na(clinic[[activePredictor]]), "Empty")

    }
  }

  #### PATIENT GROUP CREATION ####
  # Subset of Control's that remained as Control's
  cPrima <- subset(clinic,  clinic[[classVariable]] == "Control")

  # Subset of Case that changed to Control
  cPrimaPrima <- subset(clinic,  clinic[[classVariable]] == firstGroup)

  # Subset of Case's that remained as Case's
  ePrima <- subset(clinic, clinic[[classVariable]] == "Case")

  # Subset of Control that changed to Case
  ePrimaPrima <- subset(clinic,  clinic[[classVariable]] == secondGroup)

  # Total Control subset (before genetic)
  c <- rbind(cPrima, ePrimaPrima)

  # Total Case subset (before genetic)
  e <- rbind(ePrima, cPrimaPrima)

  #### RATIO ANALYSIS ####
  totalResult <- data.frame(SumCCNumerator = numeric(0), SumCCDenominator = numeric(0), SumCMNumerator = numeric(0), SumCMDenominator = numeric(0), CCCM = numeric(0), SumECNumerator = numeric(0), SumECDenominator = numeric(0), SumEMNumerator = numeric(0), SumEMDenominator = numeric(0), ECEM = numeric(0), SumENumerator = numeric(0), SumEDenominator = numeric(0), SumCNumerator = numeric(0), SumCDenominator = numeric(0), EC = numeric(0), stringsAsFactors = FALSE)

  categoricalResult <- data.frame(SumCCNumerator = numeric(0), SumCCDenominator = numeric(0), SumCMNumerator = numeric(0), SumCMDenominator = numeric(0), CCCM = numeric(0), SumECNumerator = numeric(0), SumECDenominator = numeric(0), SumEMNumerator = numeric(0), SumEMDenominator = numeric(0), ECEM = numeric(0), SumENumerator = numeric(0), SumEDenominator = numeric(0), SumCNumerator = numeric(0), SumCDenominator = numeric(0), EC = numeric(0), stringsAsFactors = FALSE)

  numericResult <- data.frame(NumPatientsCC = numeric(0), MeanCC = numeric(0), NumPatientsCM = numeric(0), MeanCM = numeric(0), CCCM = numeric(0), NumPatientsEC = numeric(0), MeanEC = numeric(0), NumPatientsEM = numeric(0), MeanEM = numeric(0), ECEM = numeric(0), NumPatientsE = numeric(0), MeanE = numeric(0), NumPatientsC = numeric(0), MeanC = numeric(0), EC = numeric(0), stringsAsFactors = FALSE)

  namesNumericRows <- character(0)
  namesCategoricRows <- character(0)
  namesAllRows <- character(0)

  # Iterating through the selected activePredictors for analysis
  for(activePredictor in activePredictors){

    # Checking if the activePredictor is numerical
    if(is.numeric(clinic[[activePredictor]])){

      # Obtaining the number of patients for the numerator and the denominator
      # Since the activePredictor is numeric, there will be empty columns
      NumPatientsCC <- nrow(ePrimaPrima)

      if(NumPatientsCC == 0){
        MeanCC <- 0
      } else {
        MeanCC <- mean(ePrimaPrima[[activePredictor]], na.rm = TRUE)
      }

      NumPatientsCM <- nrow(cPrima)

      if(NumPatientsCM == 0){
        MeanCM <- 0
      } else {
        MeanCM <- mean(cPrima[[activePredictor]], na.rm = TRUE)
      }


      NumPatientsEC <- nrow(cPrimaPrima)

      if(NumPatientsEC == 0){
        MeanEC <- 0
      } else {
        MeanEC <- mean(cPrimaPrima[[activePredictor]], na.rm = TRUE)
      }


      NumPatientsEM <- nrow(ePrima)

      if(NumPatientsEM == 0){
        MeanEM <- 0
      } else {
        MeanEM <- mean(ePrima[[activePredictor]], na.rm = TRUE)
      }


      NumPatientsC <- nrow(c)

      if(NumPatientsEM == 0){
        MeanC <- 0
      } else {
        MeanC <- mean(c[[activePredictor]], na.rm = TRUE)
      }


      NumPatientsE <- nrow(e)

      if(NumPatientsE == 0){
        MeanE <- 0
      } else {
        MeanE <- mean(e[[activePredictor]], na.rm = TRUE)
      }

      # Obtaining the desired ratios, as it is numeric, we calculate the means
      ratio1 <- MeanCC / MeanCM
      ratio1 <- replace(ratio1, is.na(ratio1), 0)
      ratio1 <- replace(ratio1, is.infinite(ratio1), 0)


      ratio2 <- MeanEC / MeanEM
      ratio2 <- replace(ratio2, is.na(ratio2), 0)
      ratio2 <- replace(ratio2, is.infinite(ratio2), 0)

      ratio3 <- MeanE / MeanC
      ratio3 <- replace(ratio3, is.na(ratio3), 0)
      ratio3 <- replace(ratio3, is.infinite(ratio3), 0)

      namesNumericRows <- c(namesNumericRows, activePredictor)

      # Creating a row for these new ratios
      newRow <- data.frame(NumPatientsCC = NumPatientsCC, MeanCC = MeanCC, NumPatientsCM = NumPatientsCM, MeanCM = MeanCM, CCCM = ratio1,
                           NumPatientsEC = NumPatientsEC, MeanEC = MeanEC, NumPatientsEM = NumPatientsEM, MeanEM = MeanEM, ECEM = ratio2,
                           NumPatientsC = NumPatientsC, MeanC = MeanC, NumPatientsE = NumPatientsE, MeanE = MeanE, CE = ratio3,
                           stringsAsFactors = FALSE)

      # Adding the new row
      numericResult <- rbind(numericResult, newRow)

      # Table for plotting
      namesAllRows <- c(namesAllRows, activePredictor)

      # Creating a row for these new ratios
      newRow <- data.frame(SumCCNumerator = NumPatientsCC, SumCCDenominator = 0, SumCMNumerator = NumPatientsCM, SumCMDenominator = 0, CCCM = ratio1,
                           SumECNumerator = NumPatientsEC, SumECDenominator = 0, SumEMNumerator = NumPatientsEM, SumEMDenominator = 0, ECEM = ratio2,
                           SumCNumerator = NumPatientsC, SumCDenominator = 0, SumENumerator = NumPatientsE, SumEDenominator = 0, CE = ratio3,
                           stringsAsFactors = FALSE)

      # Adding the new row
      totalResult <- rbind(totalResult, newRow)

    } else {

      # Predictors with a single value
      if(length(unique(clinic[[activePredictor]])) == 1){

        # Obtaining the number of patients for the numerator and the denominator
        # As there's only one value for this active predictor, there will be empty columns.
        sumCCNumerator <-  sum(ePrimaPrima[[activePredictor]] == unique(clinic[[activePredictor]]))
        sumCCDenominator <- 0

        sumCMNumerator <- sum(cPrima[[activePredictor]] == unique(clinic[[activePredictor]]))
        sumCMDenominator <- 0


        sumECNumerator <- sum(cPrimaPrima[[activePredictor]] == unique(clinic[[activePredictor]]))
        sumECDenominator <- 0

        sumEMNumerator <- sum(ePrima[[activePredictor]] == unique(clinic[[activePredictor]]))
        sumEMDenominator <- 0


        sumCNumerator<- sum(c[[activePredictor]] == unique(clinic[[activePredictor]]))
        sumCDenominator <- 0

        sumENumerator <- sum(e[[activePredictor]] == unique(clinic[[activePredictor]]))
        sumEDenominator <- 0

        # If there's only one category, the active predictor is not relevant
        # We set the ratio to 0

        ratio1 <- 0
        ratio2 <- 0
        ratio3 <- 0


        namesCategoricRows <- c(namesCategoricRows, activePredictor)

        # Create a row for these new ratios
        newRow <- data.frame(SumCCNumerator = sumCCNumerator, SumCCDenominator = sumCCDenominator, SumCMNumerator = sumCMNumerator, SumCMDenominator = sumCMDenominator, CCCM = ratio1,
                             SumECNumerator = sumECNumerator, SumECDenominator = sumECDenominator, SumEMNumerator = sumEMNumerator, SumEMDenominator = sumEMDenominator, ECEM = ratio2,
                             SumCNumerator = sumCNumerator, SumCDenominator = sumCDenominator, SumENumerator = sumENumerator, SumEDenominator = sumEDenominator, CE = ratio3,
                             stringsAsFactors = FALSE)

        # Add the new row
        categoricalResult <- rbind(categoricalResult, newRow)

        namesAllRows <- c(namesAllRows, activePredictor)

        # Adding the new row
        totalResult <- rbind(totalResult, newRow)


      # Predictors with 3 values or fewer
      } else if (length(unique(clinic[[activePredictor]])) < 4){


        # If there are 3 or fewer, we calculate all possible combinations two by two
        combinations <- combn(unique(clinic[[activePredictor]]), 2)

        # For each of these combinations
        for(i in ncol(combinations)){

          # I get the number of patients for the numerator and the denominator
          sumCCNumerator <-  sum(ePrimaPrima[[activePredictor]] == combinations[1, i])
          sumCCDenominator <- sum(ePrimaPrima[[activePredictor]] == combinations[2, i])

          sumCMNumerator <- sum(cPrima[[activePredictor]] == combinations[1, i])
          sumCMDenominator <- sum(cPrima[[activePredictor]] == combinations[2, i])

          sumECNumerator <- sum(cPrimaPrima[[activePredictor]] == combinations[1, i])
          sumECDenominator <- sum(cPrimaPrima[[activePredictor]] == combinations[2, i])

          sumEMNumerator <- sum(ePrima[[activePredictor]] == combinations[1, i])
          sumEMDenominator <- sum(ePrima[[activePredictor]] == combinations[2, i])

          sumCNumerator <- sum(c[[activePredictor]] == combinations[1, i])
          sumCDenominator <- sum(c[[activePredictor]] == combinations[2, i])

          sumENumerator <- sum(e[[activePredictor]] == combinations[1, i])
          sumEDenominator <- sum(e[[activePredictor]] == combinations[2, i])

          # Obtaining the desired ratios

          ratio1 <- (sum(ePrimaPrima[[activePredictor]] == combinations[1, i]) / sum(ePrimaPrima[[activePredictor]] == combinations[2, i])) / (sum(cPrima[[activePredictor]] == combinations[1, i]) / sum(cPrima[[activePredictor]] == combinations[2, i]))
          # If the ratio is NA/NaN (due to no patients with that activePredictor value) replace it with 0
          ratio1 <- replace(ratio1, is.na(ratio1), 0)
          # If the ratio is Infinite (due to no patients with that activePredictor value) replace it with 0
          ratio1 <- replace(ratio1, is.infinite(ratio1), 0)

          # Ratio between:
          ratio2 <- (sum(cPrimaPrima[[activePredictor]] == combinations[1, i]) / sum(cPrimaPrima[[activePredictor]] == combinations[2, i])) / (sum(ePrima[[activePredictor]] == combinations[1, i]) / sum(ePrima[[activePredictor]] == combinations[2, i]))
          # If the ratio is NA/NaN (due to no patients with that activePredictor value) replace it with 0
          ratio2 <- replace(ratio2, is.na(ratio2), 0)
          # If the ratio is Infinite (due to no patients with that activePredictor value) replace it with 0
          ratio2 <- replace(ratio2, is.infinite(ratio2), 0)

          # Ratio between:
          ratio3 <- (sum(c[[activePredictor]] == combinations[1, i]) / sum(c[[activePredictor]] == combinations[2, i])) / (sum(e[[activePredictor]] == combinations[1, i]) / sum(e[[activePredictor]] == combinations[2, i]))
          # If the ratio is NA/NaN (due to no patients with that activePredictor value) replace it with 0
          ratio3 <- replace(ratio3, is.na(ratio3), 0)
          # If the ratio is Infinite (due to no patients with that activePredictor value) replace it with 0
          ratio3 <- replace(ratio3, is.infinite(ratio3), 0)


          # We obtain a string formed by the two values separated by /
          combinationName <- paste(activePredictor, "", sep = " ")
          combinationName <- paste(combinationName, combinations[1, i], sep = "")
          combinationName <- paste(combinationName, combinations[2, i], sep = " vs ")

          namesCategoricRows <- c(namesCategoricRows, combinationName)

          # Create a row for these new ratios
          newRow <- data.frame(SumCCNumerator = sumCCNumerator, SumCCDenominator = sumCCDenominator, SumCMNumerator = sumCMNumerator, SumCMDenominator = sumCMDenominator, CCCM = ratio1,
                               SumECNumerator = sumECNumerator, SumECDenominator = sumECDenominator, SumEMNumerator = sumEMNumerator, SumEMDenominator = sumEMDenominator, ECEM = ratio2,
                               SumCNumerator = sumCNumerator, SumCDenominator = sumCDenominator, SumENumerator = sumENumerator, SumEDenominator = sumEDenominator, CE = ratio3,
                               stringsAsFactors = FALSE)

          # Add the new row
          categoricalResult <- rbind(categoricalResult, newRow)

          namesAllRows <- c(namesAllRows, combinationName)

          # Add the new row
          totalResult <- rbind(totalResult, newRow)
        }

        # Predictors with 3 values or more
      } else {

        # If not, we calculate each category against the rest
        categories <- unique(clinic[[activePredictor]])

        # For each of those categories
        for(category in categories){

          # I obtain the number of patients for the numerator and the denominator
          sumCCNumerator <-  sum(ePrimaPrima[[activePredictor]] == category)
          sumCCDenominator <- sum(ePrimaPrima[[activePredictor]] != category)

          sumCMNumerator <- sum(cPrima[[activePredictor]] == category)
          sumCMDenominator <- sum(cPrima[[activePredictor]] != category)

          sumECNumerator <- sum(cPrimaPrima[[activePredictor]] == category)
          sumECDenominator <- sum(cPrimaPrima[[activePredictor]] != category)

          sumEMNumerator <- sum(ePrima[[activePredictor]] == category)
          sumEMDenominator <- sum(ePrima[[activePredictor]] != category)

          sumCNumerator <- sum(c[[activePredictor]] == category)
          sumCDenominator <- sum(c[[activePredictor]] != category)

          sumENumerator <- sum(e[[activePredictor]] == category)
          sumEDenominator <- sum(e[[activePredictor]] != category)

          # Obtaining the desired ratios

          # Ratio between:
          ratio1 <- (sum(ePrimaPrima[[activePredictor]] == category) / sum(ePrimaPrima[[activePredictor]] != category)) / (sum(cPrima[[activePredictor]] == category) / sum(cPrima[[activePredictor]] != category))
          # If the ratio is NA/NaN (because there are no patients with that value of the active predictor), replace it with 0
          ratio1 <- replace(ratio1, is.na(ratio1), 0)
          # If the ratio is infinite (because there are no patients with that value of the active predictor), replace it with 0
          ratio1 <- replace(ratio1, is.infinite(ratio1), 0)


          # Ratio between:
          ratio2 <- (sum(cPrimaPrima[[activePredictor]] == category) / sum(cPrimaPrima[[activePredictor]] != category)) / (sum(ePrima[[activePredictor]] == category) / sum(ePrima[[activePredictor]] != category))
          # If the ratio is NA/NaN (because there are no patients with that value of the active predictor), replace it with 0
          ratio2 <- replace(ratio2, is.na(ratio2), 0)
          # If the ratio is infinite (because there are no patients with that value of the active predictor), replace it with 0
          ratio2 <- replace(ratio2, is.infinite(ratio2), 0)

          # Ratio between:
          ratio3 <- (sum(c[[activePredictor]] == category) / sum(c[[activePredictor]] != category)) / (sum(e[[activePredictor]] == category) / sum(e[[activePredictor]] != category))
          # If the ratio is NA/NaN (because there are no patients with that value of the active predictor), replace it with 0
          ratio3 <- replace(ratio3, is.na(ratio3), 0)
          # If the ratio is infinite (because there are no patients with that value of the active predictor), replace it with 0
          ratio3 <- replace(ratio3, is.infinite(ratio3), 0)

          # We obtain a string formed by the two values separated by /
          categoryName <- paste(activePredictor, "", sep = " ")
          categoryName <- paste(categoryName, category, sep = "")
          categoryName <- paste(categoryName, "Rest", sep = " vs ")

          namesCategoricRows <- c(namesCategoricRows, categoryName)


          # Create a row for these new ratios
          newRow <- data.frame(SumCCNumerator = sumCCNumerator, SumCCDenominator = sumCCDenominator, SumCMNumerator = sumCMNumerator, SumCMDenominator = sumCMDenominator, CCCM = ratio1,
                               SumECNumerator = sumECNumerator, SumECDenominator = sumECDenominator, SumEMNumerator = sumEMNumerator, SumEMDenominator = sumEMDenominator, ECEM = ratio2,
                               SumCNumerator = sumCNumerator, SumCDenominator = sumCDenominator, SumENumerator = sumENumerator, SumEDenominator = sumEDenominator, CE = ratio3,
                               stringsAsFactors = FALSE)
          # Add the new row
          categoricalResult <- rbind(categoricalResult, newRow)

          namesAllRows <- c(namesAllRows, categoryName)

          # Add the new row
          totalResult <- rbind(totalResult, newRow)
        }

      }

    }

  }

  # Establezco el nombre de las filas como el de los activePredictores
  rownames(categoricalResult) <- namesCategoricRows
  rownames(numericResult) <- namesNumericRows
  rownames(totalResult) <- namesAllRows


  # Save the ratio analysis
  solutionPathCategoric <- paste(name, "CategoricTable.tsv", sep="_")
  solutionPathNumeric <- paste(name, "NumericTable.tsv", sep="_")
  solutionPathtTotal <- paste(name, "TotalTable.tsv", sep="_")

  write.table(categoricalResult, solutionPathCategoric, row.names = T, col.names = T, sep =  '\t')
  write.table(numericResult, solutionPathNumeric, row.names = T, col.names = T, sep =  '\t')
  write.table(totalResult, solutionPathtTotal, row.names = T, col.names = T, sep =  '\t')
}
