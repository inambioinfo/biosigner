\name{getAccuracyMN}

\alias{getAccuracyMN}
\alias{getAccuracyMN,biosign-method}

\title{Accuracies of the full model and the models restricted to the signatures}

\description{
    Balanced accuracies for the full models, and the models restricted to the 'S' and 'AS' signatures
}

\usage{
\S4method{getAccuracyMN}{biosign}(object)
}

\arguments{
  \item{object}{
    An S4 object of class \code{biosign}, created by the \code{biosign} function.
  }
}

\value{
  A numeric matrix containing the balanced accuracies for the full models, and the models
  restricted to the 'S' and 'AS' signatures (predictions are obtained by
  using the resampling scheme selected with the 'bootI' and 'crossvalI'
  arguments)
}

\author{
  Philippe Rinaudo and Etienne Thevenot (CEA)
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

## individual boxplot of the selected signatures

getAccuracyMN(diaSign)

detach(diaplasma)

}