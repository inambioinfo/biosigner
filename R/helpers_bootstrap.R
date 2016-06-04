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


