\name{show.biosign}

\alias{show.biosign}
\alias{show,biosign-method}

\title{Show method for 'biosign' signature objects}

\description{
  Prints the selected features and the accuracies of the classifiers.
}

\usage{
\S4method{show}{biosign}(object)
}

\arguments{
  \item{object}{
    An S4 object of class \code{biosign}, created by the \code{biosign} function.
  }
}

\value{
  Invisible.
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

diaSign

detach(diaplasma)

}