## Project: biosigner
## Details: Contains wrappers functions
## Authors: Philippe Rinaudo and Etienne Thevenot (CEA)


######## SPECIFICATIONS: Wrapper for get.FSI function (simplify usage)
##                       Handle SVM and PLS
## INPUT: DSV: a DSV object,
##        response (character): a response name, i.e. name of a column of the observations of DSV
##        k (interger): number of folds
##        nperm (integer): number permutation
## OUTPUT: return list with the FSI, the mean rank (used for stepwise) and the score of each variable, and the evaluation of the model

fsi <- function(DSV=NULL, response=NULL,
                method=NULL,
                nboot=50,
                nperm = 1,
                alpha = 0.05, 
                fixed.rank = FALSE,
                return.full.model = FALSE){
  switch(method,
         "plsda" = {
             intern.method = "opls"
             methods.args = list(permI = 0, printL = FALSE, plotL = FALSE)
             method.x.y.name = list(x="x", y="y")
             predict.obj.newdat.name = list(object="object", newdata="newdata")
             predict.args= NULL
         },
         "randomforest" = {
             intern.method = "randomForest"
             methods.args = list(importance=TRUE,localImp=TRUE, proximity=TRUE)
             method.x.y.name = list(x="x", y="y")
             predict.obj.newdat.name = list(object="object", newdata="newdata")
             predict.args = NULL
             ## removing all variables with NA
             ind <- unique(which(is.na(DSV$dataMatrix), arr.ind = TRUE)[, "col"])
             if(length(ind)) {
                 varNamVc <- colnames(DSV$dataMatrix)[ind]
                 warning("The following variables(s) contain(s) 'NA' values and will be removed from the RandomForest models: '",
                         paste(varNamVc, collapse = "', '"),
                         "'")
                 DSV$dataMatrix <- DSV$dataMatrix[,-ind,drop=FALSE]
                 DSV$variableMetadata <- DSV$variableMetadata[-ind,,drop=FALSE]
             }
         },
         "svm" = {
             intern.method = "svm"
             methods.args = list(kernel="linear", degree=3, cost=1)
             method.x.y.name = list(x="x", y="y")
             predict.obj.newdat.name = list(object="object", newdata="newdata")
             predict.args=list(probability=FALSE)
             ## removing all variables with NA
             ind <- unique(which(is.na(DSV$dataMatrix), arr.ind = TRUE)[, "col"])
             if(length(ind)) {
                 varNamVc <- colnames(DSV$dataMatrix)[ind]
                 warning("The following variables(s) contain(s) 'NA' values and will be removed from the SVM models: '",
                         paste(varNamVc, collapse = "', '"),
                         "'")
                 DSV$dataMatrix <- DSV$dataMatrix[,-ind,drop=FALSE]
                 DSV$variableMetadata <- DSV$variableMetadata[-ind,,drop=FALSE]
             }
         },
         stop(paste0("(wrapperF.R) undefined method, must be in ('svm', 'plsda', 'randomforest'), but attempt to use ",method))
         )


  my.fsi <- get.FSI(DSV = DSV, response = response,
                    method = intern.method,
                    method.x.y.name = method.x.y.name,
                    method.args = methods.args,
                    predict.obj.newdat.name = predict.obj.newdat.name,
                    predict.args=predict.args,
                    nboot=nboot,
                    nperm = nperm,
                    alpha = alpha, 
                    fixed.rank = fixed.rank)

  if(return.full.model){
    full.model <- get.model(x = DSV$dataMatrix, y=DSV$sampleMetadata[,response],
                            method = intern.method,
                            method.x.y.name = method.x.y.name,
                            method.args = methods.args)
    my.fsi$model <- full.model
  }

  return(my.fsi)
}
