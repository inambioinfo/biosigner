#' Molecular signature discovery from omics data
#'
#' Feature selection is critical in omics data analysis to extract restricted
#' and meaningful molecular signatures from complex and high-dimension data, and
#' to build robust classifiers. This package implements a new method to assess
#' the relevance of the variables for the prediction performances of the
#' classifier. The approach can be run in parallel with the PLS-DA, Random
#' Forest, and SVM binary classifiers. The signatures and the corresponding
#' 'restricted' models are returned, enabling future predictions on new
#' datasets. A Galaxy implementation of the package is available within the
#' Workflow4metabolomics.org online infrastructure for computational
#' metabolomics.
#'
#' @import methods e1071 randomForest ropls
#' @importFrom Biobase ExpressionSet exprs pData
#' @importFrom grDevices dev.new dev.off pdf png
#' @importFrom graphics abline arrows axis box boxplot image layout mtext par rect title
#' @importFrom stats median var
#' @importFrom utils head tail
#'
#' @name biosigner-package
#' @aliases biosigner-package biosigner
#' @docType package
#' @author Philippe Rinaudo <phd.rinaudo@@gmail.com> and Etienne Thevenot
#' <etienne.thevenot@@cea.fr>.
#'
#' Maintainer: Philippe Rinaudo <phd.rinaudo@@gmail.com>
#' @keywords package
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
NULL


