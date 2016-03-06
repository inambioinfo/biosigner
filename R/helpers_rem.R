## Project: biosigner
## Details: core functions to generate and manipulate models from DSV objects
##          Mainly used REM object for Rank, evaluation and model.
##          REM also contains score and feature
## Authors: Philippe Rinaudo and Etienne Thevenot (CEA)


## INPUT: x: a matrix of numerics (row: observations, column: features), used to train the model
##        y: a vector of responses (number of responses must match the number of observations/row of x), used to train the model
##        method: a R function, with at least 2 inputs: a numerical matrix and a vector of responses
##        method.x.y.name: arguments name for the input matrix and the input vector of responses of the method (x,y by default)
##        method.args: a list of arguments for the R function "method"
get.model <- function(x=NULL, y=NULL,
                      method=NULL,
                      method.x.y.name = list(x="x", y="y"),
                      method.args = NULL){
  ## generate the full arguments for the model training
  #   generate the names of the args (to match the x and y)
  args.name.full <- c(method.x.y.name$x,
                      method.x.y.name$y,
                      names(method.args))
  #   generate the args list with the arg names just created
  args.full <- c(list(x,y),method.args)
  names(args.full) <- args.name.full
  ## generate the model from training set
  model <- do.call(method, args.full)
  return(model)
}

## INPUT: model: object containing a model
##        x.test: a matrix of numerics (row: observations, column: features), used to test the model
##        predict.obj.newdat.name: list of arguments names for the input model and the input newdata of responses of the prediction method (object,newdata by default)
##        predict.args: a list of arguments for the prediction "method" (except model and newdata)
get.prediction <- function(model=NULL,
                           x.test=NULL,
                           predict.obj.newdat.name = list(object="object", newdata="newdata"),
                           predict.args=NULL){

  ## generate the prediction
  ## generate the full arguments for the model prediction
  args.name.full <- c(predict.obj.newdat.name$object,
                      predict.obj.newdat.name$newdata,
                      names(predict.args))
  #   generate the args list with the arg names just created
  args.full <- c(list(model,x.test),predict.args)
  names(args.full) <- args.name.full

  ## generate the prediction
  prediction <- do.call(predict, args.full)

  return(prediction)
}

## INPUT: prediction: a vector of prediction
##        y.test: a vector of responses
## OUTPUT: evaluation: evaluation of the prediction: compute the balanced classification rate
##                      for regression: R square
##                      for classification: balanced accuracy
get.evaluation <- function(prediction = NULL,
                           y.test = NULL){
    if(class(y.test) != "numeric"){
      y.test <- as.character(y.test)
      if(class(prediction) == "numeric")
        stop("predicted outpout is numeric and original output is not numeric (factor or character)")
      prediction <- as.character(prediction)
      class.y.test <- unique(y.test)
      if(length(class.y.test) != 2)
        stop("Trying to estimate balanced classification rate with more or less than two classes")
      rate <- 0
      for(i in class.y.test){
        ind <- which(y.test == i)
        rate <- rate + (sum(prediction[ind] == y.test[ind]) / length(ind))
      }
      rate <- rate / length(class.y.test)
      return(rate)
    }
    else{
        SSE <- sum( (y.test-prediction)^2 )
        SST <- sum( (y.test-mean(y.test))^2 )
        evaluation <- 1 - SSE/SST
        return(evaluation)
    }
}

## INPUT: model: object containing a model
##        x.train: the training set on which the model has been trained
##        x.test: a matrix of numerics (row: observations, column: features), used to test the model
##        y.test: response to predict with x.test as data
##        evaluation: evaluation of the model generate with x.test and y.test
##        predict.obj.newdat.name: list of arguments names for the input model and the input newdata of responses of the prediction method (object,newdata by default)
##        predict.args: a list of arguments for the prediction "method" (except model and newdata)
get.var.score <- function(model = NULL,
                          x.train = NULL,
                          x.test = NULL,
                          y.test = NULL,
                          evaluation = NULL,
                          predict.obj.newdat.name = list(object="object", newdata="newdata"),
                          predict.args=NULL){
  ## manual handling for svm:
  if(class(model)=="svm"){
    score <- (t(model$coefs) %*% x.train[model$index,])^2
  }
  else if(class(model)=="opls"){
    #initialised the score
    score <- rep(0, ncol(x.train))
    ## get the variable with a vip:
    ind.vip <- which(colnames(x.train) %in% names(model$vipVn))
    score[ind.vip] <- model$vipVn
  }
  else if(class(model)=="randomForest"){
    score <- model$importance[,1]
  }
  else{
    ## generate the score for each features, TODO: generalize this given a score criteria in input (make a function), for example make a variable importance function
    score <- c()
    for(i in 1:ncol(x.train)){
      x.test.perm <- x.test
      x.test.perm[,i] <- sample(x.test[,i])
      prediction.perm <- get.prediction(model,x.test.perm,predict.obj.newdat.name,predict.args)
      score <- c(score, evaluation - get.evaluation(prediction.perm,y.test))
    }
  }
  names(score) <- colnames(x.train)
  return(score)
}

######## SPECIFICATIONS: train a model from a training data set and a training method (with args)
########                 then test the model on a test set and compute the accuracy
########                 then evaluate the feature importance by feature permutation
## INPUT: x: a matrix of numerics (row: observations, column: features), used to train the model
##        y: a vector of responses (number of responses must match the number of observations/row of x), used to train the model
##        method: a R function, with at least 2 inputs: a numercial matrix and a vector of responses
##        method.x.y.name: list of arguments names for the input matrix and the input vector of responses of the method (x,y by default)
##        method.args: a list of arguments for the R function "method" (except data and response)
##        x.test: a matrix of numerics (row: observations, column: features), used to test the model
##        y.test: a vector of responses (number of responses must match the number of observations/row of x.test), used to test the model
##        predict.obj.newdat.name: list of arguments names for the input model and the input newdata of responses of the prediction method (object,newdata by default)
##        predict.args: a list of arguments for the prediction "method" (except model and newdata)
##        TODO: add an argument to choose how to evaluate the model (currently accuracy only)
## OUTPUT: a REM, containing the rank of the features, their scores, the evaluation of the model and the model.
get.REM <- function(x=NULL, y=NULL,
                    method=NULL,
                    method.x.y.name = list(x="x", y="y"),
                    method.args = NULL,
                    x.test=NULL, y.test=NULL,
                    predict.obj.newdat.name = list(object="object", newdata="newdata"),
                    predict.args=NULL){
  # removing cte variables
  n <- ncol(x)
  var.names <- colnames(x)
  var.cte <- which(apply(x, 2, var) <= .Machine["double.eps"])
  if(length(var.cte) > 0)
    x <- x[, -var.cte, drop=FALSE]

  model <- get.model(x,y,method,method.x.y.name,method.args)

  if(length(var.cte) > 0)
    x.test <- x.test[, -var.cte, drop=FALSE]
  
  prediction <- get.prediction(model,x.test,predict.obj.newdat.name,predict.args)
  ## generate the evaluation
  evaluation <- get.evaluation(prediction, y.test)

  ## generate the score
  score <- get.var.score(model = model,
                         x.train = x,
                         x.test = x.test, y.test = y.test,
                         evaluation = evaluation,
                         predict.obj.newdat.name = predict.obj.newdat.name, predict.args = predict.args)
  score.full <- rep(0, n)
  if(length(var.cte > 0))
    score.full[-var.cte] <- score
  else
    score.full <- score
  ## generate the rank
  rank <- rank(-score.full, ties.method="max")
  names(rank) <- var.names

  # generate the REM
  REM <- list(rank=rank, evaluation=evaluation, model=model, score=score.full, features=var.names, var.cte=var.cte)
  return(REM)
}

######## SPECIFICATIONS: generate the mean evaluation, score and correspoing rank from a list of REMs (generate from a same DSV)
## INPUT: REMs: a list of REM objects,
## OUTPUT: a list containing the mean evaluation, score and correspoing rank and features
get.REMs.summary <- function(REMs=NULL){
  nb.rem <- length(REMs)
  ## get the cumulative evaluation and scores
  # evaluation
  cumul.evaluation <- 0
  nb.feat <- length(REMs[[1]]$score)
  score <- matrix(nrow=length(REMs), ncol=nb.feat)

  for(i in 1:length(REMs)){
    cumul.evaluation <- cumul.evaluation + REMs[[i]]$evaluation
    score[i,] <- rank(-REMs[[i]]$score)
  }
  evaluation <- cumul.evaluation / nb.rem
  score <- apply(score, 2, median)
  rank <- rank(score, ties.method="max")

  return(list(rank=rank, evaluation=evaluation, score=score, feature=REMs[[1]]$feature))
}

######## SPECIFICATIONS: generate a list of REM from bootstrap  or cross validationsampling, see get.REM for more details
## INPUT: DSV: a DSV object,
##        response: the response name (column of the sampleMetadata)
##        nboot: number of bootstrap (for bootstrap)
##        k: number of folds (for kfold)
##        see get.REM for the other arguments
## OUTPUT: a list of "nboot" REM
get.REMs <- function(DSV=NULL, response=NULL,
                    method=NULL,
                    method.x.y.name = list(x="x", y="y"),
                    method.args = NULL,
                    predict.obj.newdat.name = list(object="object", newdata="newdata"),
                    predict.args=NULL,
                    nboot=0,
                    fixed.rank = FALSE){

  ## generate the all REM
  REM.list <- list()
  if(nboot!=0){
    for(i in 1:nboot){
      # generate the bootstrap
      DSV <- extract.bootstrap(DSV, boot.sunder = TRUE, response = response)
      # generate the corresponding REM
      REM <- get.REM(x = get.training.profile.boot(DSV), y = get.training.response.boot(DSV,response),
                      method = method,
                      method.x.y.name = method.x.y.name,
                      method.args = method.args,
                      x.test = get.test.profile.boot(DSV), y.test = get.test.response.boot(DSV,response),
                      predict.obj.newdat.name = predict.obj.newdat.name,
                      predict.args = predict.args)
      REM$ind.test <- get.ind.test.profile.boot(DSV)
      # add to the heap
      REM.list <- c(REM.list,list(REM))
    }
  }

  ## generate general evaluation, score...
  REM.summary <- get.REMs.summary(REM.list)
  if(fixed.rank){
    fixed.model <- get.model(x = DSV$dataMatrix, y = DSV$sampleMetadata$response,
                                  method = method,
                                  method.x.y.name = method.x.y.name,
                                  method.args = method.args)
    REM.summary$rank <- get.var.score(model = fixed.model,
                                      x.train =  DSV$dataMatrix,
                                      x.test =  DSV$dataMatrix, y.test = DSV$sampleMetadata$response,
                                      evaluation = REM.summary$evaluation,
                                      predict.obj.newdat.name = predict.obj.newdat.name,
                                      predict.args = predict.args)
  }

  REMs <- list(DSV = DSV, REM=REM.list, rank=REM.summary$rank, evaluation=REM.summary$evaluation, score=REM.summary$score, feature=REM.summary$feature)

  return(REMs)
}

get.REMs.FSI <- function(REMs=NULL, DSV=NULL, response=NULL,
                         predict.obj.newdat.name = list(object="object", newdata="newdata"),
                         predict.args=NULL,
                         nperm = 0,
                         alpha = 0.05){
  #research the FSI the closest to alpha by dichotomy
  ind.order <- order(REMs$rank)
  bound.min <- 1
  bound.max <- length(ind.order)
  stop.criteria <- FALSE
  i <- bound.min
  fsi.bound.min <- 1
  while( (bound.min < bound.max) | (bound.min == bound.max & length(ind.order) == 1) ){
    evaluation <- rep(0,nperm*length(REMs$REM))
    # generate the permutated data nperm times
    ind <- which(REMs$rank >= REMs$rank[ind.order[i]])
    for( j in 1:length(REMs$REM)){  
      for( s in 1:nperm ){
        # permute variables
        profile.perm <- apply(DSV$dataMatrix[,ind,drop=FALSE], 2, sample)
        profile <- DSV$dataMatrix
        profile[,ind] <- profile.perm
        # generate prediction model
        # get test profile
        if(length(REMs$REM[[j]]$var.cte) > 0){
          profile.rem <- profile[REMs$REM[[j]]$ind.test, -REMs$REM[[j]]$var.cte, drop=FALSE]
        }
        else{
          profile.rem <- profile[REMs$REM[[j]]$ind.test, , drop=FALSE]
        }
        # predict
        prediction <- get.prediction(REMs$REM[[j]]$model,profile.rem,predict.obj.newdat.name,predict.args)
        # generate evaluation
        y.test <- DSV$sampleMetadata[,response][REMs$REM[[j]]$ind.test]
        evaluation[(s-1)*length(REMs$REM)+j] <- REMs$REM[[j]]$evaluation - get.evaluation(prediction = prediction, y.test = y.test)
      }
    }
    fsi.current <- sum( evaluation <= 0 )/(length(REMs$REM)*nperm)
    if(fsi.current <= alpha){
      bound.min <- i
      fsi.bound.min <- fsi.current
    }
    else
      bound.max <- i-1
    if(length(ind.order)==1)
      bound.max <- bound.min - 1
    i <- ceiling((bound.min+bound.max) / 2)
  }
  FSI <- rep(-1, length(REMs$rank))
  names(FSI) <- rownames(DSV$variableMetadata)
  if(fsi.bound.min <= alpha)
    FSI[ind.order[1:bound.min]] <- 0
  return(FSI)
}

get.FSI <- function(DSV=NULL, response=NULL,
                     method=NULL,
                     method.x.y.name = list(x="x", y="y"),
                     method.args = NULL,
                     predict.obj.newdat.name = list(object="object", newdata="newdata"),
                     predict.args=NULL,
                     nboot=0,
                     nperm = 0,
                     alpha = 0, 
                     fixed.rank = FALSE){
  REMs <- get.REMs(DSV = DSV, response = response,
                   method=method,
                   method.x.y.name = method.x.y.name,
                   method.args = method.args,
                   predict.obj.newdat.name = predict.obj.newdat.name,
                   predict.args=predict.args,
                   nboot=nboot,
                   fixed.rank = fixed.rank)
  if(nperm > 0){
    FSI <- get.REMs.FSI(REMs=REMs, DSV=REMs$DSV, response=response,
                      predict.obj.newdat.name = predict.obj.newdat.name,
                      predict.args = predict.args,
                      nperm = nperm,
                      alpha = alpha)
    REMs$FSI <- FSI
  }
  return(REMs)
}






