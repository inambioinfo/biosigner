## Package: biosigner
## Details: functions to generate and manipulate bootstrap from DSV objects
## Authors: Philippe Rinaudo and Etienne Thevenot (CEA)


## INPUT: DSV: a DSV object,
##        boot.size (interger): the size of the number of obervations to sample (default is set to the number of observations of DSV)
##        boot.sunder (boolean): controle the proportion in the boostrapping
##        response (character): a response name, i.e. name of a column of the observations of DSV
## OUTPUT: return a the DSV input with a supplementary dataframe (bootDF) contaning all needed information to extract train and test sub DSV.
extract.bootstrap <- function(DSV=NULL, boot.size=0, boot.sunder=FALSE, response=NULL){
  # Initialization of the bootDF
  bootDF <- data.frame(row.names=row.names(DSV$sampleMetadata))
  bootDF$Freq <- 0
  if(boot.size==0)
    boot.size <- nrow(DSV$sampleMetadata)

  # Bootstrapping of the observations
  boot.sampling <- c()
  if(boot.sunder & !is.null(response)){
    # get the responses
    # for all responses get the relative number of sample by bootstrapping
    boot.sampling <- c()
    res.all <- unique(DSV$sampleMetadata[[response]])
    for(res in res.all){
      prop.res <- sum( DSV$sampleMetadata[[response]] == res )
      # Remove 0.2 samples to ensure 0.2 OOB samples
      obs.res <- DSV$sampleMetadata[which(DSV$sampleMetadata[,response]==res),,drop=FALSE]
      ind.0.8 <- sort(sample(x = 1:nrow(obs.res), size = round(0.8*nrow(obs.res)), replace = FALSE)) 
      boot.sampling <- c(boot.sampling,table(sample(row.names(obs.res[ind.0.8,,drop=FALSE]), size=prop.res, replace=TRUE)))
    }
  }
  else
    boot.sampling <- table(sample(row.names(DSV$sampleMetadata), size=boot.size, replace=TRUE))

  # Updating the bootDF
  indices <- match(names(boot.sampling), rownames(bootDF))

  bootDF$Freq[indices] <- boot.sampling
  # Set Train or not in the Train column
  bootDF$Train <- bootDF$Freq != 0

  # Checking for constant variables

  # Add the bootDF to the input DSV object and return the DSV object
  DSV$bootDF <- bootDF
  return(DSV)
}

## INPUT: DSV: a DSV object,
## OUTPUT: a profile matrix, generated from the bootstrap dataframe of DSV and corresponding to the observations tagged as TRAIN=TRUE with their frequencies
get.training.profile.boot <- function(DSV=NULL){
  if(is.null(DSV$bootDF))
    return(DSV$dataMatrix)
  # get profile with only training observations
  train.sub.profile <- DSV$dataMatrix[which(DSV$bootDF$Train==TRUE),,drop=FALSE]
  train.boot <- DSV$bootDF[which(DSV$bootDF$Train==TRUE),,drop=FALSE]

  boot.size <- sum(DSV$bootDF$Freq)

  # Duplicating profile sample relatively to the bootstrap freq
  train.profile <- matrix(rep(train.sub.profile, rep(train.boot$Freq,ncol(train.sub.profile))),boot.size, ncol(train.sub.profile))
  colnames(train.profile) <- colnames(DSV$dataMatrix)
  return(train.profile)
}

## INPUT: DSV: a DSV object,
##        response (character): a response name, i.e. name of a column of the observations of DSV
## OUTPUT: a factor response, generated from the bootstrap dataframe of DSV and corresponding to the observations tagged as TRAIN=TRUE with their frequencies
get.training.response.boot <- function(DSV=NULL, response=NULL){
  if(is.null(DSV$bootDF))
    return(DSV$sampleMetadata[,response])
  # get only training observations
  train.sub.obs<- DSV$sampleMetadata[which(DSV$bootDF$Train==TRUE),,drop=FALSE][,response]
  train.boot <- DSV$bootDF[which(DSV$bootDF$Train==TRUE),,drop=FALSE]

  # Duplicating sample response relatively to the bootstrap freq
  train.response <- factor(rep(train.sub.obs, train.boot$Freq))
  return(train.response)
}

## INPUT: DSV: a DSV object,
## OUTPUT: a profile matrix, generated from the bootstrap dataframe of DSV and corresponding to the observations tagged as TRAIN=FALSE
get.test.profile.boot <- function(DSV=NULL){
  if(is.null(DSV$bootDF))
    return(DSV$dataMatrix)
  # get profile with only test observations
  test.profile <- DSV$dataMatrix[which(DSV$bootDF$Train==FALSE),,drop=FALSE]

  return(test.profile)
}

## INPUT: DSV: a DSV object,
##        response (character): a response name, i.e. name of a column of the observations of DSV
## OUTPUT: a factor response, generated from the bootstrap dataframe of DSV and corresponding to the observations tagged as TRAIN=FALSE
get.test.response.boot <- function(DSV=NULL, response=NULL){
  if(is.null(DSV$bootDF))
    return(DSV$sampleMetadata[,response])
  # get only test observations
  test.response <- factor(DSV$sampleMetadata[which(DSV$bootDF$Train==FALSE),,drop=FALSE][,response])
  return(test.response)
}

get.ind.test.profile.boot <- function(DSV=NULL){
  return(which(DSV$bootDF$Train==FALSE))
}

