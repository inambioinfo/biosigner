\name{biosign}
\docType{methods}

\alias{biosign}
\alias{biosign,data.frame-method}
\alias{biosign,matrix-method}

\title{
  Builds the molecular signature.
}

\description{
  Main function of the 'biosigner' package. For each of the available
  classifiers (PLS-DA, Random Forest, and SVM), the significant features
  are selected and the corresponding models are built.
}

\usage{
\S4method{biosign}{data.frame}(x, ...)

\S4method{biosign}{matrix}(x,
y,
methodVc = c("all", "plsda", "randomforest", "svm")[1],
bootI = 50,
pvalN = 0.05,

permI = 1,
fixRankL = FALSE,

printL = TRUE,
plotL = TRUE,

.sinkC = NULL,
...)
}

\arguments{
  \item{x}{
    Numerical data frame or matrix (observations x variables); NAs are allowed for
    PLS-DA but for SVM, samples with NA will be removed
  }
  \item{y}{
    Factor: Two-level factor corresponding to the class labels
  }
  \item{methodVc}{
    Character vector: Either one or all of the following classifiers:
    Partial Least Squares Discriminant Analysis ('plsda'), or Random Forest ('randomforest'), or Support Vector Machine ('svm')
  }
  \item{bootI}{
    Integer: Number of bootstaps for resampling
  }
  \item{pvalN}{
    Numeric: To speed up the selection, only variables which
    significantly improve the model up to two times this threshold (to
    take into account potential fluctuations) are computed
  }
  \item{permI}{
    Integer: Random permutation are used to assess the significance of each new variable included into the model (forward selection)
  }
  \item{fixRankL}{
    Logical: Should the initial ranking be computed with the full model only,
    or as the median of the ranks from the models built on the sampled dataset?
  }
  \item{printL}{
    Logical: Should informations regarding the data set and the model be
    printed? [default = TRUE]
  }
  \item{plotL}{
    Logical: Should the 'summary' plot be displayed? [default = TRUE]
  }
  \item{.sinkC}{
    Character: Name of the file for R output diversion [default = NULL:
    no diversion]; Diversion of messages is required for the integration
    into Galaxy
  }
  \item{\dots}{
    Currently not used.
  }
}

\value{
  An S4 object of class 'biosign' containing the following slots: 1) 'methodVc' character vector: selected classifier(s)
  ('plsda', 'randomforest', and/or 'svm'), 2) 'accuracyMN' numeric
  matrix: balanced accuracies for the full models, and the models
  restricted to the 'S' and 'AS' signatures (predictions are obtained by
  using the resampling scheme selected with the 'bootI' and 'crossvalI'
  arguments), 3) 'tierMC' character matrix: contains the tier ('S',
  'A', 'B', 'C', 'D', or 'E') of each feature for each classifier
  (features with tier 'S' have been found significant in all
  backward selections; features with tier 'A' have been found
  significant in all but the last selection, and so on), 4)
  modelLs list: selected classifier(s) trained on the subset restricted to
  the 'S' features, 5) signatureLs list: 'S' signatures for each
  classifier; and 6) 'AS' list: 'AS' signatures and corresponding
  trained classifiers, in addition to the dataset restricted to tiers
  'S' and 'A' ('xMN') and the labels ('yFc')
}

\author{
Philippe Rinaudo and Etienne Thevenot (CEA)
}

\seealso{\code{\link{predict.biosign}}, \code{\link{plot.biosign}}}


\examples{

## loading the diaplasma dataset

data(diaplasma)
attach(diaplasma)

## restricting to a smaller dataset for this example

featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
dataMatrix <- dataMatrix[, featureSelVl]
variableMetadata <- variableMetadata[featureSelVl, ]

## signature selection for all 3 classifiers
## a bootI = 5 number of bootstraps is used for this example
## we recommend to keep the default bootI = 50 value for your analyzes

set.seed(123)
diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)

detach(diaplasma)

}

