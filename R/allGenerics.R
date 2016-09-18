#' Builds the molecular signature.
#'
#' Main function of the 'biosigner' package. For each of the available
#' classifiers (PLS-DA, Random Forest, and SVM), the significant features are
#' selected and the corresponding models are built.
#'
#' @name biosign
#' @rdname biosign
#' @aliases biosign biosign,data.frame-method biosign,matrix-method
#' @docType methods
#' @param x Numerical data frame or matrix (observations x variables), or
#' ExpressionSet object with non empty assayData and phenoData; NAs are allowed
#' for PLS-DA but for SVM, samples with NA will be removed
#' @param y Two-level factor corresponding to the class labels, or a character
#' indicating the name of the column of the pData to be used, when x is an
#' ExpressionSet object
#' @param methodVc Character vector: Either one or all of the following
#' classifiers: Partial Least Squares Discriminant Analysis ('plsda'), or
#' Random Forest ('randomforest'), or Support Vector Machine ('svm')
#' @param bootI Integer: Number of bootstaps for resampling
#' @param pvalN Numeric: To speed up the selection, only variables which
#' significantly improve the model up to two times this threshold (to take into
#' account potential fluctuations) are computed
#' @param permI Integer: Random permutation are used to assess the significance
#' of each new variable included into the model (forward selection)
#' @param fixRankL Logical: Should the initial ranking be computed with the
#' full model only, or as the median of the ranks from the models built on the
#' sampled dataset?
#' @param printL Logical: Should informations regarding the data set and the
#' model be printed? [default = TRUE]
#' @param plotL Logical: Should the 'summary' plot be displayed? [default =
#' TRUE]
#' @param .sinkC Character: Name of the file for R output diversion [default =
#' NULL: no diversion]; Diversion of messages is required for the integration
#' into Galaxy
#' @param ... Currently not used.
#' @return An S4 object of class 'biosign' containing the following slots: 1)
#' 'methodVc' character vector: selected classifier(s) ('plsda',
#' 'randomforest', and/or 'svm'), 2) 'accuracyMN' numeric matrix: balanced
#' accuracies for the full models, and the models restricted to the 'S' and
#' 'AS' signatures (predictions are obtained by using the resampling scheme
#' selected with the 'bootI' and 'crossvalI' arguments), 3) 'tierMC' character
#' matrix: contains the tier ('S', 'A', 'B', 'C', 'D', or 'E') of each feature
#' for each classifier (features with tier 'S' have been found significant in
#' all backward selections; features with tier 'A' have been found significant
#' in all but the last selection, and so on), 4) modelLs list: selected
#' classifier(s) trained on the subset restricted to the 'S' features, 5)
#' signatureLs list: 'S' signatures for each classifier; and 6) 'AS' list: 'AS'
#' signatures and corresponding trained classifiers, in addition to the dataset
#' restricted to tiers 'S' and 'A' ('xMN') and the labels ('yFc')
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @seealso \code{\link{predict.biosign}}, \code{\link{plot.biosign}}
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
#' @export
setGeneric("biosign", function(x, ...) standardGeneric("biosign"))

#' Accuracies of the full model and the models restricted to the signatures
#'
#' Balanced accuracies for the full models, and the models restricted to the
#' 'S' and 'AS' signatures
#'
#' @aliases getAccuracyMN getAccuracyMN,biosign-method
#' @param object An S4 object of class \code{biosign}, created by the
#' \code{biosign} function.
#' @param ... Currently not used.
#' @return A numeric matrix containing the balanced accuracies for the full
#' models, and the models restricted to the 'S' and 'AS' signatures
#' (predictions are obtained by using the resampling scheme selected with the
#' 'bootI' and 'crossvalI' arguments)
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
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
#' ## individual boxplot of the selected signatures
#'
#' getAccuracyMN(diaSign)
#'
#' detach(diaplasma)
#'
#' @rdname getAccuracyMN
#' @export
setGeneric("getAccuracyMN", function(object, ...) {standardGeneric("getAccuracyMN")})

#' Signatures selected by the models
#'
#' List of 'S' (or 'S' and 'A') signatures for each classifier
#'
#' @aliases getSignatureLs getSignatureLs,biosign-method
#' @param object An S4 object of class \code{biosign}, created by the
#' \code{biosign} function.
#' @param tierC Character: defines whether signatures from the 'S' tier only
#' (default) or the ('S' and 'A') tiers should be returned
#' @param ... Currently not used.
#' @return List of 'S' (or 'S' and 'A') signatures for each classifier
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
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
#' ## individual boxplot of the selected signatures
#'
#' getSignatureLs(diaSign)
#'
#' detach(diaplasma)
#'
#' @rdname getSignatureLs
#' @export
setGeneric("getSignatureLs", function(object, tierC = c("S", "AS")[1], ...) {standardGeneric("getSignatureLs")})
