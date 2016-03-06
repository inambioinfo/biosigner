\name{getSignatureLs}

\alias{getSignatureLs}
\alias{getSignatureLs,biosign-method}

\title{Signatures selected by the models}

\description{
    List of 'S' (or 'S' and 'A') signatures for each classifier
}

\usage{
\S4method{getSignatureLs}{biosign}(object, tierC = c("S", "AS")[1])
}

\arguments{
  \item{object}{
    An S4 object of class \code{biosign}, created by the \code{biosign} function.
  }
  \item{tierC}{
    Character: defines whether signatures from the 'S' tier only
    (default) or the ('S' and 'A') tiers should be returned
  }
}

\value{
  List of 'S' (or 'S' and 'A') signatures for each classifier
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

getSignatureLs(diaSign)

detach(diaplasma)

}