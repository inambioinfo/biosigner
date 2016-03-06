\name{biosign-class}
\docType{class}

\alias{biosign-class}

\title{Class "biosign"}

\description{The biosigner object class}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("biosign", ...)}
  or by calling the \code{biosign} function
}

\section{Slots}{
  \describe{
    \item{\code{methodVc}:}{character vector: selected classifier(s) ('plsda', 'randomforest', or 'svm')}
    \item{\code{accuracyMN}:}{numeric matrix: balanced accuracies for
  the full models, and the models restricted to the 'S' and 'AS' signatures}
    \item{\code{tierMC}:}{character matrix: contains the tier ('S', 'A',
      'B', 'C', 'D', or 'E') of each feature for each classifier}
    \item{\code{yFc}:}{factor with two levels: response factor}
    \item{\code{modelLs}:}{list: selected classifier(s) trained on the subset restricted to the 'S' features}
    \item{\code{signatureLs}:}{list: 'S' signatures for each classifier}
    \item{\code{xSubMN}:}{matrix: dataset restricted to the 'S' tier}    
    \item{\code{AS}:}{list: 'AS' signatures and corresponding trained classifiers, in addition to the dataset restricted to tiers 'S' and 'A' ('xMN') and the labels ('yFc')}
  }
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(object = "biosign")}: plots the tiers or the boxplots of the biomarkers selected by all the classifiers}
    \item{predict}{\code{signature(object = "biosign")}: predicts values
  based upon the classifiers  trained on the 'S' signature}
  }
}

\author{
Philippe Rinaudo and Etienne Thevenot (CEA)
}

\seealso{
  \code{\link{biosign}},
}

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

