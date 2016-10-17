#' Class "biosign"
#'
#' The biosigner object class
#'
#' @name biosign-class
#' @rdname biosign-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("biosign", ...)} or by calling the \code{biosign} function
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @seealso \code{\link{biosign}}
#'
#' @slot methodVc character vector: selected classifier(s) ('plsda', 'randomforest', or 'svm')
#' @slot accuracyMN numeric matrix: balanced accuracies for the full models, and the models restricted to the 'S' and 'AS' signatures
#' @slot tierMC character matrix: contains the tier ('S', 'A', 'B', 'C', 'D', or 'E') of each feature for each classifier
#' @slot yFc factor with two levels: response factor
#' @slot modelLs list: selected classifier(s) trained on the subset restricted to the 'S' features
#' @slot signatureLs list: 'S' signatures for each classifier
#' @slot xSubMN matrix: dataset restricted to the 'S' tier
#' @slot AS list: 'AS' signatures and corresponding trained classifiers, in addition to the dataset restricted to tiers 'S' and 'A' ('xMN') and the labels ('yFc')
#'
#' @examples
#'
#' ## loading the diaplasma dataset
#'
#' data(diaplasma)
#' attach(diaplasma)
#'
#' ## restricting to a smaller dataset for this example
#'
#' featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
#' dataMatrix <- dataMatrix[, featureSelVl]
#' variableMetadata <- variableMetadata[featureSelVl, ]
#'
#' ## signature selection for all 3 classifiers
#' ## a bootI = 5 number of bootstraps is used for this example
#' ## we recommend to keep the default bootI = 50 value for your analyzes
#'
#' set.seed(123)
#' diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)
#'
#' detach(diaplasma)
#'
#' @exportClass biosign
setClass("biosign",
         slots = c(methodVc = "character",
             accuracyMN = "matrix",
             tierMC = "matrix",
             yFc = "factor",
             modelLs = "list",
             signatureLs = "list",
             xSubMN = "matrix",
             AS = "list"))


