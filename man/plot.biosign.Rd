\name{plot.biosign}

\alias{plot.biosign}
\alias{plot,biosign-method}

\title{Plot method for 'biosign' signature objects}

\description{
    Displays classifier tiers or individual boxplots from selected features
}

\usage{
\S4method{plot}{biosign}(x,
y,
tierMaxC = "S",
typeC = c("tier", "boxplot")[1],
file.pdfC = NULL,
.sinkC = NULL,
...)
}

\arguments{
  \item{x}{
    An S4 object of class \code{biosign}, created by the \code{biosign} function.
  }
  \item{y}{
    Currently not used.
  }
  \item{tierMaxC}{
    Character: Maximum level of tiers to display: Either 'S' and 'A',
    (for boxplot), or also 'B', 'C', 'D', and 'E' (for tiers) by decreasing number of selections
  }
  \item{typeC}{
    Character: Plot type; either 'tier' [default] displaying the
    comparison of signatures up to the selected 'tierMaxC' or 'boxplot'
    showing the individual boxplots of the features selected by all the
    classifiers
  }
  \item{file.pdfC}{
    Character: Figure filename ending with '.pdf'; default is NULL (no saving; displaying instead)
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
  A plot is created on the current graphics device.
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

plot(diaSign, typeC = "boxplot")

detach(diaplasma)

}