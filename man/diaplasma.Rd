\name{diaplasma}

\alias{diaplasma}

\docType{data}

\title{
  Analysis of plasma from diabetic patients by LC-HRMS
}

\description{
  Plasma samples from 69 diabetic patients were analyzed by
  reversed-phase liquid chromatography coupled to high-resolution mass
  spectrometry (Orbitrap Exactive) in the negative ionization mode. The
  raw data were pre-processed with XCMS and CAMERA (5,501 features),
  corrected for signal drift, log10 transformed, and annotated with an
  in-house spectral database. The patient's age, body mass index, and
  diabetic type are recorded. These three clinical covariates are strongly
  associated, most of the type II patients being older and with a higher bmi
  than the type I individuals.
}

\usage{
data(diaplasma)
}

\format{
  A list with the following elements: 1) dataMatrix: a 69 samples x
  5,501 features matrix of numeric type corresponding to the intensity
  profiles (values have been log10-transformed); 2) sampleMetadata: a 69 x 3
  data frame, with the patients' diabetic type ('type', factor), age
  ('age', numeric), and body mass index ('bmi', numeric), and 3) variableMetadata: a 5,501 x 8 data
  frame, with the median m/z ('mzmed', numeric) and the median retention
  time in seconds ('rtmed', numeric) from XCMS, the 'isotopes'
  (character), 'adduct' (character) and 'pcgroups' (numeric) annotations
  from CAMERA, and the names of the m/z and
  RT matching compounds from an in-house database of pure spectra from commercial
  metabolites ('spiDb', character).
}

\value{
  List containing the 'dataMatrix' matrix (numeric) of data (samples as
  rows, variables as columns), the 'sampleMetadata' data frame of sample
  metadata, and the variableMetadata data frame of variable metadata. Row
  names of 'dataMatrix' and 'sampleMetadata' are identical. Column names of
  'dataMatrix' are identical to row names of 'variableMetadata'. For
  details see the 'Format' section above.
}

\references{
  Rinaudo, P., Junot C., and Thevenot E.A. (submitted). biosigner: A new method for the discovery of restricted and stable molecular signatures from omics data.
}

\keyword{datasets}
