## Note: Hierarchy of helper function calls:
##
##1..biosign                           methods-biosign_class
##   2..getTierF                       helpers_significance
##      3..getBootSignificanceF        helpers_significance
##         4..getBootModelF            helpers_model
##            5..getBootExtract        helpers_bootstrap
##            5..getModelAccuRankF     helpers_model
##               6..getBootTrainxF     helpers_bootstrap
##               6..getBootTrainyF     helpers_bootstrap
##               6..getBootTestxF      helpers_bootstrap
##               6..getBootTestyF      helpers_bootstrap
##               6..getModelF          helpers_model
##               6..getPredictionF     helpers_model
##               6..getAccuracyF       helpers_model
##               6..getImportanceF     helpers_model
##            5..getBootTestIndF       helpers_bootstrap
##            5..getBootSummaryF       helpers_model
##         4..getSignificanceF         helpers_significance
##            5..getPredictionF        helpers_model
##            5..getAccuracyF          helpers_model
##   2..getBootSignificanceF           helpers_significance
##      3..getModelF                   helpers_model


##-------------------##
## helpers_bootstrap ##
##-------------------##


## Package: biosigner
## Details: functions to generate and manipulate bootstrap from dataLs objects
## Authors: Philippe Rinaudo and Etienne Thevenot (CEA)


## INPUT: dataLs: a dataLs object,
##        bootI (interger): the size of the number of obervations to sample
##                          (default is set to the number of observations of dataLs)
##        bootSubSampL (boolean): Should only 80% of the samples be used for bootstraping
##                                (ensuring 20% of out of the bag samples)
##        respC (character): a response name, i.e. name of a column of the observations
##                           of dataLs
## OUTPUT: returns the dataLs input with a supplementary dataframe (bootDF) contaning
##         all needed information to extract train and test sub dataLs.
bootExtractF <- function(dataLs = NULL,
                         bootI = 0,
                         bootSubSampL = FALSE,
                         respC = NULL){

    ## Initialization of bootDF
    bootDF <- data.frame(row.names = row.names(dataLs$sampleMetadata))
    bootDF$freqVn <- numeric(nrow(bootDF))

    if(bootI == 0)
        bootI <- nrow(dataLs$sampleMetadata)

    ## Bootstrapping of the observations
    bootFreqVi <- numeric()
    if(bootSubSampL && !is.null(respC)){
        ## gets the responses
        ## for all responses gets the relative number of samples by bootstrapping
        classVc <- unique(dataLs$sampleMetadata[[respC]])

        ## classI <- min(table(dataLs$sampleMetadata[[respC]]))

        for(classC in classVc){

            classI <- sum(dataLs$sampleMetadata[[respC]] == classC)

            ## Removing 0.2 samples to ensure 0.2 OOB samples
            sampleMetaClassDF <- dataLs$sampleMetadata[which(dataLs$sampleMetadata[, respC] == classC), , drop = FALSE]
            bootClassVi <- sort(sample(1:nrow(sampleMetaClassDF),
                                       size = round(0.8 * nrow(sampleMetaClassDF)),
                                       replace = FALSE))
            bootFreqVi <- c(bootFreqVi,
                            table(sample(row.names(sampleMetaClassDF[bootClassVi, , drop = FALSE]),
                                         size = classI,
                                         replace = TRUE)))
        }
    } else
        bootFreqVi <- table(sample(row.names(dataLs$sampleMetadata),
                                   size = bootI,
                                   replace = TRUE))

    ## Updating the bootDF
    bootVi <- match(names(bootFreqVi), rownames(bootDF))
    bootDF$freqVn[bootVi] <- bootFreqVi

    ## Set TRUE or FALSE in the trainVl column
    bootDF$trainVl <- bootDF$freqVn != 0

    ## Checking for constant variables

    ## Adding the bootDF to the dataLs object
    dataLs$bootDF <- bootDF

    return(dataLs)
}

## INPUT: dataLs: a dataLs object,
## OUTPUT: a profile matrix, generated from the bootstrap dataframe of dataLs and corresponding to the observations tagged as TRAIN=TRUE with their frequencies
getBootTrainxF <- function(dataLs = NULL){

  if(is.null(dataLs$bootDF))
    return(dataLs$dataMatrix)

  ## get profile with only training observations
  xTrainMN <- dataLs$dataMatrix[which(dataLs$bootDF$trainVl == TRUE), , drop = FALSE]
  bootTrainDF <- dataLs$bootDF[which(dataLs$bootDF$trainVl == TRUE), , drop = FALSE]

  bootI <- sum(dataLs$bootDF$freqVn)

  ## Duplicating profile sample relatively to the bootstrap freq
  xTrainMN <- matrix(rep(xTrainMN,
                         rep(bootTrainDF$freqVn,
                             ncol(xTrainMN))),
                     bootI,
                     ncol(xTrainMN))
  colnames(xTrainMN) <- colnames(dataLs$dataMatrix)
  return(xTrainMN)
}

## INPUT: dataLs: a dataLs object,
##        respC (character): a response name, i.e. name of a column of the observations of dataLs
## OUTPUT: a factor response, generated from the bootstrap dataframe of dataLs and corresponding to the observations tagged as TRAIN=TRUE with their frequencies
getBootTrainyF <- function(dataLs=NULL, respC=NULL){
  if(is.null(dataLs$bootDF))
    return(dataLs$sampleMetadata[, respC])
  # get only training observations
  yTrainFc <- dataLs$sampleMetadata[which(dataLs$bootDF$trainVl == TRUE), , drop = FALSE][, respC]
  bootTrainDF <- dataLs$bootDF[which(dataLs$bootDF$trainVl == TRUE), , drop = FALSE]

  # Duplicating sample response relatively to the bootstrap freq
  yTrainFc <- factor(rep(yTrainFc, bootTrainDF$freqVn))
  return(yTrainFc)
}

## INPUT: dataLs: a dataLs object,
## OUTPUT: a profile matrix, generated from the bootstrap dataframe of dataLs and corresponding to the observations tagged as TRAIN=FALSE
getBootTestxF <- function(dataLs=NULL){

  if(is.null(dataLs$bootDF))
    return(dataLs$dataMatrix)
  # get profile with only test observations
  xTestMN <- dataLs$dataMatrix[which(dataLs$bootDF$trainVl == FALSE), , drop=FALSE]

  return(xTestMN)

}

## INPUT: dataLs: a dataLs object,
##        respC (character): a response name, i.e. name of a column of the observations of dataLs
## OUTPUT: a factor response, generated from the bootstrap dataframe of dataLs and corresponding to the observations tagged as TRAIN=FALSE
getBootTestyF <- function(dataLs=NULL, respC=NULL){

  if(is.null(dataLs$bootDF))
    return(dataLs$sampleMetadata[, respC])
  # get only test observations
  yTestFc <- factor(dataLs$sampleMetadata[which(dataLs$bootDF$trainVl==FALSE), , drop=FALSE][, respC])

  return(yTestFc)

}

getBootTestIndF <- function(dataLs = NULL)
  return(which(dataLs$bootDF$trainVl == FALSE))


##---------------##
## helpers_model ##
##---------------##


## Project: biosigner
## Details: Helper functions to build models on dataLs subsets generated by bootstraping
##          Computes the accuracy of the models, and the rank of the features
## Authors: Philippe Rinaudo and Etienne Thevenot (CEA)


## INPUT:  predTestFc: a factor of prediction
##         yTestFc: a factor of responses
## OUTPUT: accuracyN: accuracy of the prediction:
##                      for regression: R square (not used)
##                      for classification: balanced accuracy
getAccuracyF <- function(predTestFc = NULL,
                         yTestFc = NULL){

    if(class(yTestFc) != "numeric"){

      yTestVc <- as.character(yTestFc)
      yClassVc <- unique(yTestVc)
      if(length(yClassVc) != 2)
        stop("Computing accuracy for 2-class discrimination only")

      if(class(predTestFc) == "numeric")
        stop("Predicted output is numeric and original output is not numeric (factor or character)")
      predTestVc <- as.character(predTestFc)

      accuracyN <- 0
      for(yClassC in yClassVc){
        yClassVi <- which(yTestVc == yClassC)
        accuracyN <- accuracyN + (sum(predTestVc[yClassVi] == yTestVc[yClassVi]) / length(yClassVi))
      }
      accuracyN <- accuracyN / length(yClassVc)
      return(accuracyN)
    }
    ## else{ ## not used currently (for future regression applications)
    ##     SSE <- sum( (yTestFc - predTestFc)^2 )
    ##     SST <- sum( (yTestFc - mean(yTestFc))^2 )
    ##     accuracyN <- 1 - SSE/SST
    ##     return(accuracyN)
    ## }

} ## getAccuracyF


## SPECIFICATIONS: generates a list of models from bootstrap,
##                 see getModelAccuRankF for more details
## INPUT: dataLs: data and metadata,
##        respC: the response name (column of the sampleMetadata)
##        bootI: number of bootstraps
##        see getModelAccuRankF for the other arguments
## OUTPUT: a list of bootI models with accuracy and variable ranking
getBootModelF <- function(dataLs = NULL,
                          respC = NULL,
                          methC = NULL,
                          methNamLs = list(x = "x", y = "y"),
                          methArgLs = NULL,
                          predNamLs = list(object = "object", newdata = "newdata"),
                          predArgLs = NULL,
                          bootI = 0,
                          fixRankL = FALSE){

  ## generate the all the modelAccuRank generated by bootstrap
  modelAccuRankLs <- list()

  if(bootI != 0){

    for(i in 1:bootI){

      # generate the bootstrap
      dataLs <- bootExtractF(dataLs,
                             bootSubSampL = TRUE,
                             respC = respC)

      # generate the corresponding modelAccuRank
      modelAccuRank <- getModelAccuRankF(getBootTrainxF(dataLs),
                                         getBootTrainyF(dataLs, respC),
                                         methC = methC,
                                         methNamLs = methNamLs,
                                         methArgLs = methArgLs,
                                         xTestMN = getBootTestxF(dataLs),
                                         yTestFc = getBootTestyF(dataLs, respC),
                                         predNamLs = predNamLs,
                                         predArgLs = predArgLs)

      modelAccuRank$ind.test <- getBootTestIndF(dataLs)
      ## add to the heap
      modelAccuRankLs <- c(modelAccuRankLs, list(modelAccuRank))
    }
  }

  ## generate general evaluation, score...
  summaryLs <- getBootSummaryF(modelAccuRankLs)
  if(fixRankL){
      fixed.model <- getModelF(dataLs$dataMatrix,
                               dataLs$sampleMetadata$respC,
                               methC = methC,
                               methNamLs = methNamLs,
                               methArgLs = methArgLs)
      summaryLs$rankVn <- getImportanceF(fixed.model,
                                         dataLs$dataMatrix)
                                         ## dataLs$dataMatrix,
                                         ## dataLs$sampleMetadata$respC,
                                         ## summaryLs$accuracyN,
                                         ## predNamLs,
                                         ## predArgLs)
  }

  bootModelLs <- list(dataLs = dataLs,
                      modelAccuRankLs = modelAccuRankLs,
                      rankVn = summaryLs$rankVn,
                      accuracyN = summaryLs$accuracyN,
                      varImpVn = summaryLs$varImpVn,
                      varNamVc = summaryLs$varNamVc)

  return(bootModelLs)

} ## getBootModelF


## SPECIFICATIONS: generates the mean accuracy, median importance
##                 and corresponding rank from a list of modelAccuRank
##                 (generated from a same dataLs)
## INPUT: bootModelLs: a list of models with accuracy and variable rankings,
##                 generated by bootstrap (getBootModelF)
## OUTPUT: a list containing the mean evaluation, importance
##         and corresponding rank and feature names
getBootSummaryF <- function(bootModelLs = NULL){

    bootI <- length(bootModelLs)
    ## get the cumulative evaluation and scores evaluation
    accuCumN <- 0
    varI <- length(bootModelLs[[1]]$varImpVn)
    varImpMN <- matrix(nrow = bootI, ncol = varI)

    for(i in 1:bootI){
        accuCumN <- accuCumN + bootModelLs[[i]]$accuracyN
        varImpMN[i, ] <- rank(-bootModelLs[[i]]$varImpVn)
    }
    accuracyN <- accuCumN / bootI
    varImpVn <- apply(varImpMN, 2, median)
    rankVn <- rank(varImpVn, ties.method="max")

    return(list(rankVn = rankVn,
                accuracyN = accuracyN,
                varImpVn = varImpVn,
                varNamVc = bootModelLs[[1]]$varNamVc))
}


## INPUT: model: object containing a model
##        xTrainMN: the training set on which the model has been trained
##        xTestMN: a matrix of numerics (row: observations, column: features), used to test the model
##        yTestFc: response to predict with xTestMN as data
##        accuracyN: evaluation of the model generated with xTestMN and yTestFc
##        predNamLs: list of arguments names for the input model and the input newdata of responses of the prediction method (object, newdata by default)
##        predArgLs: a list of arguments for the prediction "methC" (except model and newdata)
## OUTPUT: numerical vector of feature importance
getImportanceF <- function(model = NULL,
                           xTrainMN = NULL
                           ## xTestMN = NULL,
                           ## yTestFc = NULL,
                           ## accuracyN = NULL,
                           ## predNamLs = list(object = "object",
                           ##     newdata = "newdata"),
                           ## predArgLs=NULL
                           ){

    if(class(model) == "svm"){ ## weights (w)

        varImpVn <- (t(model[["coefs"]]) %*% xTrainMN[model[["index"]],])^2

    }
    else if(class(model) == "opls"){ ## VIP

        ## initialize
        varImpVn <- rep(0, ncol(xTrainMN))
        ## get the variable with a VIP:
        vipVn <- getVipVn(model)
        vipVi <- which(colnames(xTrainMN) %in% names(vipVn))
        varImpVn[vipVi] <- vipVn

    }
    else if(class(model) == "randomForest"){

        varImpVn <- model[["importance"]][, 1]

    }
    ## else{ ## not currently used
    ##     ## generic metric to measure variable importance
    ##     ## TODO: generalize this given a score criteria in input (make a function), for example make a variable importance function
        ## varImpVn <- numeric(ncol(xTrainMN))
        ## for(j in 1:length(varImpVn)){
        ##     xTestMN.perm <- xTestMN
        ##     xTestMN.perm[, j] <- sample(xTestMN[, j])
        ##     predTestFc.perm <- getPredictionF(model,
        ##                                       xTestMN.perm,
        ##                                       predNamLs,
        ##                                       predArgLs)
        ##     varImpVn[j] <- accuracyN - getAccuracyF(predTestFc.perm,
        ##                                             yTestFc)
        ## }
    ## }
    names(varImpVn) <- colnames(xTrainMN)

    return(varImpVn)

} ## getImportanceF


## INPUT: xMN: a numerical matrix of predictors (row: observations, column: features),
##             used to train the model
##        yFc: a response factor (number of responses must match the number of
##             observations/row of xMN), used to train the model
##        methC: name of the classifier, i.e. an R function, with at least 2 inputs:
##               a numerical matrix and a response factor
##        methNamLs: arguments name for the input matrix and the input factor of
##                   responses of the method (x, y by default)
##        methArgLs: a list of arguments for the 'methC' R classifier
getModelF <- function(xMN = NULL,
                      yFc = NULL,
                      methC = NULL,
                      methNamLs = list(x = "x", y = "y"),
                      methArgLs = NULL){
    ## generating the full arguments for the model training
    ## generating the names of the args (to match the x and y)
    argNamFulVc <- c(methNamLs[["x"]],
                     methNamLs[["y"]],
                     names(methArgLs))

    ## generating the args list with the arg names just created
    argFulLs <- c(list(xMN, yFc), methArgLs)
    names(argFulLs) <- argNamFulVc
    ## generates the model from training set
    if(methC == "opls") {
        model <- try(do.call(methC, argFulLs), silent = TRUE)
        if(inherits(model, "try-error") &&
           substr(unclass(attr(model, "condition"))$message, 1, 85) == "No model was built because the first predictive component was already not significant") {
            argFulLs <- c(argFulLs, list(predI = 1))
            model <- do.call(methC, argFulLs)
        }
    } else
        model <- do.call(methC, argFulLs)

    return(model)
}


## SPECIFICATIONS: trains a model from a training data set and a training method (with args)
##                 then tests the model on a test set and compute the accuracy
##                 then evaluates the feature importance by feature permutation
## INPUT: xMN: a matrix of numerics (row: observations, column: features),
##             used to train the model
##        yFc: a response factor (number of responses must match the number of
##             observations/row of xMN), used to train the model
##        methC: a R function, with at least 2 inputs: a numercial matrix and a
##               vector of responses
##        methNamLs: list of arguments names for the input matrix and the input
##                   vector of responses of the method (x,y by default)
##        methArgLs: a list of arguments for the R function "methC" (except data and
##                   response)
##        xTestMN: a matrix of numerics (row: observations, column: features),
##                 used to test the model
##        yTestFc: a vector of responses (number of responses must match the number
##                 of observations/row of xTestMN), used to test the model
##        predNamLs: list of arguments names for the input model and the input newdata
##                   of responses of the prediction method (object,newdata by default)
##        predArgLs: a list of arguments for the prediction "methC"
##                   (except model and newdata)
##        TODO: add an argument to choose how to evaluate the model
##              (currently accuracy only)
## OUTPUT: a list containing the rank of the features, their importance,
##         the accuracy of the model and the model itself.
getModelAccuRankF <- function(xMN = NULL,
                              yFc = NULL,
                              methC = NULL,
                              methNamLs = list(x = "x",
                                  y = "y"),
                              methArgLs = NULL,
                              xTestMN = NULL,
                              yTestFc = NULL,
                              predNamLs = list(object = "object",
                                  newdata = "newdata"),
                              predArgLs = NULL){
  # removing constant variables
  varI <- ncol(xMN)
  varNamVc <- colnames(xMN)
  varCstVi <- which(apply(xMN, 2, var) <= .Machine["double.eps"])
  if(length(varCstVi) > 0)
    xMN <- xMN[, -varCstVi, drop = FALSE]

  model <- getModelF(xMN = xMN,
                     yFc = yFc,
                     methC = methC,
                     methNamLs = methNamLs,
                     methArgLs = methArgLs)

  if(length(varCstVi) > 0)
    xTestMN <- xTestMN[, -varCstVi, drop = FALSE]

  predTestFc <- getPredictionF(model = model,
                               xTestMN = xTestMN,
                               predNamLs = predNamLs,
                               predArgLs = predArgLs)
  ## accuracy
  accuracyN <- getAccuracyF(predTestFc = predTestFc,
                            yTestFc = yTestFc)

  ## feature importance

  varImpVn <- getImportanceF(model = model,
                             xTrainMN = xMN)
                             ## xTestMN = xTestMN,
                             ## yTestFc = yTestFc,
                             ## accuracyN = accuracyN,
                             ## predNamLs = predNamLs,
                             ## predArgLs = predArgLs)
  varImpFulVn <- rep(0, varI)
  if(length(varCstVi > 0))
    varImpFulVn[-varCstVi] <- varImpVn
  else
    varImpFulVn <- varImpVn
  ## generate the rank
  rankVn <- rank(-varImpFulVn, ties.method = "max")
  names(rankVn) <- varNamVc

  # generate the REM
  modelAccuRank <- list(rankVn = rankVn,
                        accuracyN = accuracyN,
                        model = model,
                        varImpVn = varImpFulVn,
                        varNamVc = varNamVc,
                        varCstVi = varCstVi)
  return(modelAccuRank)

} ## getModelAccuRankF


## INPUT: model: object containing a model (classifier)
##        xTestMN: numerical matrix corresponding to the test subset (row: observations, column: features)
##        predNamLs: list of arguments names for the input model and the input newdata of responses of the prediction method (object, newdata by default)
##        predArgLs: a list of arguments for the prediction "methC" (in addition to the object model and newdata)
## OUTPUT: factor of predictions
getPredictionF <- function(model = NULL,
                           xTestMN = NULL,
                           predNamLs = list(object = "object",
                               newdata = "newdata"),
                           predArgLs = NULL){

    ## generates the prediction
    ## generates the full arguments for the model prediction
    argNamFulVc <- c(predNamLs$object,
                     predNamLs$newdata,
                     names(predArgLs))
    ## generates the args list with the arg names just created
    argFulLs <- c(list(model,xTestMN),predArgLs)
    names(argFulLs) <- argNamFulVc

    ## generates the prediction
    predTestFc <- do.call(predict, argFulLs)

    return(predTestFc)

} ## getPredictionF


##----------------------##
## helpers_significance ##
##----------------------##


## Project: biosigner
## Details: Helper functions to evaluate the significance of the features
##          by using the models generated by bootstrap
##          and to iterate the process to compute the tiers
## Authors: Philippe Rinaudo and Etienne Thevenot (CEA)


## SPECIFICATIONS: For a given method (PLS-DA, Random Forest, and SVM)
##                 at a given iteration step (i.e. tier; dataset restricted to
##                                            a selected number of features)
##                 computes the bootstrap models, the aggregated accuracy
##                 and the variable significance
## INPUT: dataLs: data and metadata
##        respC (character): a response name, i.e. name of a column of the
##                           observations of dataLs
##        permI (integer): number permutations
## OUTPUT: list with the FSI, the mean rank (used for stepwise) and the significance of each variable, and the accuracy of the model
getBootSignificanceF <- function(dataLs = NULL,
                                 respC = NULL,
                                 methC = NULL,
                                 bootI = 50,
                                 permI = 1,
                                 pvalN = 0.05,
                                 fixRankL = FALSE,
                                 fullModelL = FALSE){

    switch(methC,
           "plsda" = {
               methC = "opls"
               methArgLs = list(permI = 0, printL = FALSE, plotL = FALSE)
               methNamLs = list(x="x", y="y")
               predNamLs = list(object="object", newdata="newdata")
               predArgLs= NULL
           },
           "randomforest" = {
               methC = "randomForest"
               methArgLs = list(importance=TRUE,localImp=TRUE, proximity=TRUE)
               methNamLs = list(x="x", y="y")
               predNamLs = list(object="object", newdata="newdata")
               predArgLs = NULL
               ## removing all variables with NA
               ind <- unique(which(is.na(dataLs$dataMatrix), arr.ind = TRUE)[, "col"])
               if(length(ind)) {
                   varNamVc <- colnames(dataLs$dataMatrix)[ind]
                   warning("The following variables(s) contain(s) 'NA' values and will be removed from the RandomForest models: '",
                           paste(varNamVc, collapse = "', '"),
                           "'")
                   dataLs$dataMatrix <- dataLs$dataMatrix[,-ind,drop=FALSE]
                   dataLs$variableMetadata <- dataLs$variableMetadata[-ind,,drop=FALSE]
               }
           },
           "svm" = {
               methC = "svm"
               methArgLs = list(kernel="linear", cost=1)
               ## methArgLs = list(kernel="linear", cost=10)
               ## methArgLs = list(kernel="radial", cost=10)
               methNamLs = list(x="x", y="y")
               predNamLs = list(object="object", newdata="newdata")
               predArgLs = list(probability = FALSE)
               ## removing all variables with NA
               ind <- unique(which(is.na(dataLs$dataMatrix), arr.ind = TRUE)[, "col"])
               if(length(ind)) {
                   varNamVc <- colnames(dataLs$dataMatrix)[ind]
                   warning("The following variables(s) contain(s) 'NA' values and will be removed from the SVM models: '",
                           paste(varNamVc, collapse = "', '"),
                           "'")
                   dataLs$dataMatrix <- dataLs$dataMatrix[,-ind,drop=FALSE]
                   dataLs$variableMetadata <- dataLs$variableMetadata[-ind,,drop=FALSE]
               }
           },
           stop(paste0("Undefined method, must be in ('svm', 'plsda', 'randomforest'), but attempt to use ", methC))
           )

    bootModelLs <- getBootModelF(dataLs = dataLs,
                                 respC = respC,
                                 methC = methC,
                                 methNamLs = methNamLs,
                                 methArgLs = methArgLs,
                                 predNamLs = predNamLs,
                                 predArgLs = predArgLs,
                                 bootI = bootI,
                                 fixRankL = fixRankL)

    if(permI > 0){

        signifVn <- getSignificanceF(bootModelLs = bootModelLs,
                                     dataLs = bootModelLs$dataLs,
                                     respC = respC,
                                     predNamLs = predNamLs,
                                     predArgLs = predArgLs,
                                     permI = permI,
                                     pvalN = pvalN)

        bootModelLs$signifVn <- signifVn

    }

    if(fullModelL){
        full.model <- getModelF(dataLs$dataMatrix,
                                dataLs$sampleMetadata[, respC],
                                methC = methC,
                                methNamLs = methNamLs,
                                methArgLs = methArgLs)
        bootModelLs$model <- full.model
    }

    return(bootModelLs)

} ## getBootSignificanceF


## Looks for the significant feature closest to pvalN by dichotomy
## returns a numeric vector indicating for each variable if it is relevant (0) or not (-1)
getSignificanceF <- function(bootModelLs = NULL,
                             dataLs = NULL,
                             respC = NULL,
                             predNamLs = list(object = "object",
                                 newdata = "newdata"),
                             predArgLs = NULL,
                             permI = 0,
                             pvalN = 0.05){

    bootI <- length(bootModelLs$modelAccuRankLs)

    rankOrdVi <- order(bootModelLs$rankVn)

    boundMinI <- 1
    boundMaxI <- length(rankOrdVi)

    ## stopL <- FALSE

    i <- boundMinI

    boundMinSignifN <- 1

    while((boundMinI < boundMaxI) |
          (boundMinI == boundMaxI & length(rankOrdVi) == 1)){

        accuracyVn <- rep(0, permI * bootI)

        ## generate the permutated data permI times
        highRankVi <- which(bootModelLs$rankVn >= bootModelLs$rankVn[rankOrdVi[i]])

        for(bI in 1:bootI){ ## number of bootstraps

            for(pI in 1:permI){ ## number of permutations
                ## permute variables

                xTrainPermMN <- dataLs$dataMatrix
                xTrainPermMN[, highRankVi] <- apply(dataLs$dataMatrix[, highRankVi, drop=FALSE], 2, sample)

                ## generate prediction model
                ## get test profile
                if(length(bootModelLs$modelAccuRankLs[[bI]]$varCstVi) > 0){
                    xTestPermMN <- xTrainPermMN[bootModelLs$modelAccuRankLs[[bI]]$ind.test,
                                                -bootModelLs$modelAccuRankLs[[bI]]$varCstVi, drop = FALSE]
                } else
                    xTestPermMN <- xTrainPermMN[bootModelLs$modelAccuRankLs[[bI]]$ind.test, , drop=FALSE]

                ## predict
                predTestPermFc <- getPredictionF(bootModelLs$modelAccuRankLs[[bI]]$model,
                                                 xTestPermMN,
                                                 predNamLs = predNamLs,
                                                 predArgLs = predArgLs)
                ## generate evaluation
                yTestFc <- dataLs$sampleMetadata[, respC][bootModelLs$modelAccuRankLs[[bI]]$ind.test]
                accuracyVn[(pI - 1) * bootI + bI] <- bootModelLs$modelAccuRankLs[[bI]]$accuracyN - getAccuracyF(predTestPermFc,
                                                                                                                yTestFc = yTestFc)
                ## negative if the model prediction on the randomized test subset
                ## performs better than on the true test subset
            }
        }

        signifN <- sum(accuracyVn <= 0, na.rm = TRUE) / (bootI * permI)
        if(signifN <= pvalN){
            ## i.e. predictions on the randomized test subset rarely outperforms the true test subset
            ## i.e. some of the permutated features are relevant
            ## i.e. the searched interval is shifted to the features of highest ranks
            boundMinI <- i
            boundMinSignifN <- signifN
        }
        else
            boundMaxI <- i - 1

        if(length(rankOrdVi) == 1)
            boundMaxI <- boundMinI - 1

        i <- ceiling((boundMinI + boundMaxI) / 2)

    }

    signifVn <- rep(-1, length(bootModelLs$rankVn))
    names(signifVn) <- rownames(dataLs$variableMetadata)

    if(boundMinSignifN <= pvalN)
        signifVn[rankOrdVi[1:boundMinI]] <- 0

    return(signifVn)

} ## getSignificanceF


getTierF <- function(datasetLs,
                     methodC,
                     bootI,
                     permI,
                     pvalN,
                     fixRankL) {

    fsiResLs <- list()

    tierVn <- numeric(ncol(datasetLs[["dataMatrix"]]))
    names(tierVn) <- colnames(datasetLs[["dataMatrix"]])

    stopIterL <- FALSE

    fsiResLs <- getBootSignificanceF(datasetLs,
                                     "respC",
                                     methodC,
                                     bootI = bootI,
                                     permI = permI,
                                     pvalN = pvalN,
                                     fixRankL = fixRankL)
    ## Begin iteration step
    ## extract significant variable in the last recursive step
    varSelVc <- names(which(fsiResLs$signifVn <= pvalN &
                            fsiResLs$signifVn >= 0))
    if(length(varSelVc) < 1){
        varSelVc <- fsiResLs$varNamVc[which(fsiResLs$rankVn %in% 1:(nrow(datasetLs$variableMetadata)/2))]
        varSelVi <- which(rownames(datasetLs$variableMetadata) %in% varSelVc)
        tierVn[varSelVi] <- tierVn[varSelVi] - 1
    }
    stopIterL <- FALSE
    nb.rec <- 0
    ## fsi.rec <- rep(-1, nrow(datasetLs$variableMetadata))
    ## Next line is not actually used
    ## permI.rec <- permI * ceiling(nrow(datasetLs$variableMetadata) / length(varSelVc))
    while(length(varSelVc) >= 1 && stopIterL == FALSE){
        ## while the significant variables of the current step is not the same as the ones of the previous step (stop.rec)
        ## extract the data with only the significant variable of the previous step
        datasetLs.rec <- datasetLs
        varSelVi <- which(rownames(datasetLs$variableMetadata) %in% varSelVc)
        ## increase the tiers of the previous significant variables
        tierVn[varSelVi] <- tierVn[varSelVi] + 1
        ## extraction
        datasetLs.rec$dataMatrix <- datasetLs.rec$dataMatrix[, varSelVi, drop = FALSE]
        datasetLs.rec$variableMetadata <- datasetLs.rec$variableMetadata[varSelVi, , drop = FALSE]
        ## computing the fsi for current step
        fsi.rec <- getBootSignificanceF(datasetLs.rec,
                                        "respC",
                                        methodC,
                                        bootI = bootI,
                                        permI = permI,
                                        pvalN = pvalN,
                                        fixRankL = fixRankL)
        ## checking if the number of significant variables has changed, else stop the rec
        if(length(varSelVc) != length(names(which(fsi.rec$signifVn <= pvalN &
                     fsi.rec$signifVn >= 0))) ){
            ## check if we need to iterate on the best half part (ie no significant variable found).
            if(length(varSelVc) > 1 & length(names(which(fsi.rec$signifVn <= pvalN & fsi.rec$signifVn >= 0))) < 1 ){
                varSelVc <- fsi.rec$varNamVc[which(fsi.rec$rankVn %in% 1:(nrow(datasetLs.rec$variableMetadata)/2))]
                varSelVi <- which(rownames(datasetLs$variableMetadata) %in% varSelVc)
                tierVn[varSelVi] <- tierVn[varSelVi] - 1
            }
            else{
                varSelVc <- names(which(fsi.rec$signifVn <= pvalN & fsi.rec$signifVn >= 0))
            }
            ## permI.rec <- permI.rec * ceiling(nrow(datasetLs.rec$variableMetadata) / length(varSelVc))
        }
        else{
            stopIterL <- TRUE
        }
    }
    ## update the fsiResLs TODO: ACTUALLY NOT USED ANYMORE !! CHECK IT
    varSelVi <- which(rownames(datasetLs$variableMetadata) %in% varSelVc)
    fsiResLs$signifVn[varSelVi] <- fsi.rec$signifVn
    fsiResLs$signifVn[-varSelVi] <- -1

    return(list(tierVn = tierVn,
                accuracyN = fsiResLs$accuracyN,
                stopIterL = stopIterL))

} ## getTierF
