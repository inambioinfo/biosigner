\name{predict.biosign}

\alias{predict.biosign}
\alias{predict,biosign-method}

\title{
    Predict method for 'biosign' signature objects
}

\description{
   This function predicts values based upon \code{biosign} classifiers
   trained on the 'S' signature
}

\usage{
\S4method{predict}{biosign}(object, newdata, tierMaxC = "S", ...)
}

\arguments{
  \item{object}{
    An S4 object of class \code{biosign}, created by \code{biosign} function.
  }
  \item{newdata}{
    Either a data frame or a matrix, containing numeric columns only,
    with column names identical to the 'x' used
    for model training with 'biosign'.
  }
  \item{tierMaxC}{
    Character: Maximum level of tiers to display: Either 'S'or 'A'.
  }    
  \item{\dots}{
    Currently not used.
  }
}

\value{
  Data frame with the predictions for each classifier as factor columns.
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

## training the classifiers
## a bootI = 5 number of bootstraps is used for this example
## we recommend to keep the default bootI = 50 value for your analyzes

set.seed(123)
diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)

## fitted values (for the subsets restricted to the 'S' signatures)
sFitDF <- predict(diaSign)

## confusion tables
print(lapply(sFitDF, function(predFc) table(actual = sampleMetadata[,
"type"], predicted = predFc)))

## balanced accuracies
sapply(sFitDF, function(predFc) { conf <- table(sampleMetadata[,
"type"], predFc)
conf <- sweep(conf, 1, rowSums(conf), "/")
mean(diag(conf))
})
## note that these values are slightly different from the accuracies
## returned by biosign because the latter are computed by using the
## resampling scheme selected by the bootI or crossvalI arguments
getAccuracyMN(diaSign)["S", ]

detach(diaplasma)
}