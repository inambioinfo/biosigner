\name{biosigner-package}

\alias{biosigner-package}

\alias{biosigner}

\docType{package}

\title{
  Molecular signature discovery from omics data
}

\description{
  Molecular signature discovery from omics data
}

\details{
  \tabular{ll}{
    Package: \tab biosigner\cr
    Type: \tab Package\cr
    Version: \tab 0.99.9\cr
    Date: \tab 2016-03-04\cr
    License: \tab CeCILL\cr
    LazyLoad: \tab yes\cr
  }
}

\author{
  Philippe Rinaudo <phd.rinaudo@gmail.com> and Etienne Thevenot <etienne.thevenot@cea.fr>.

  Maintainer: Philippe Rinaudo <phd.rinaudo@gmail.com>

}

\keyword{ package }

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
