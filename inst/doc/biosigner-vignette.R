## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=6, fig.path='figures/')

## ----loading, message = FALSE--------------------------------------------
library(biosigner)

## ----diaplasma-----------------------------------------------------------
data(diaplasma)

## ----diaplasma_strF------------------------------------------------------
attach(diaplasma)
library(ropls)
strF(dataMatrix)
strF(sampleMetadata)
strF(variableMetadata)

## ----diaplasma_plot, eval=FALSE------------------------------------------
#  with(sampleMetadata,
#  plot(age, bmi, cex = 1.5, col = ifelse(type == "T1", "blue", "red"), pch = 16))
#  legend("topleft", cex = 1.5, legend = paste0("T", 1:2),
#  text.col = c("blue", "red"))

## ----diaplasma_plot_fig, echo = FALSE------------------------------------
with(sampleMetadata,
plot(age, bmi, cex = 1.5, col = ifelse(type == "T1", "blue", "red"), pch = 16))
legend("topleft", cex = 1.5, legend = paste0("T", 1:2),
text.col = c("blue", "red"))

## ----select--------------------------------------------------------------
featureSelVl <- variableMetadata[, "mzmed"] >= 450 &
variableMetadata[, "mzmed"] < 500
sum(featureSelVl)
dataMatrix <- dataMatrix[, featureSelVl]
variableMetadata <- variableMetadata[featureSelVl, ]

## ----biosign_false, eval = FALSE-----------------------------------------
#  set.seed(123)
#  diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)
#  set.seed(NULL)

## ----biosign, echo = FALSE-----------------------------------------------
set.seed(123)
diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5,
plotL = FALSE)
plot(diaSign, file.pdfC = "figures/diasign.png")
set.seed(NULL)

## ----boxplot_false, eval = FALSE-----------------------------------------
#  plot(diaSign, typeC = "boxplot")

## ----boxplot, echo = FALSE-----------------------------------------------
plot(diaSign, typeC = "boxplot", file.pdfC = "figures/diasign_box.png")

## ----signature-----------------------------------------------------------
variableMetadata[getSignatureLs(diaSign)[["complete"]], ]

## ----train---------------------------------------------------------------
trainVi <- 1:floor(0.8 * nrow(dataMatrix))
testVi <- setdiff(1:nrow(dataMatrix), trainVi)

## ----biosign_train_false, eval = FALSE-----------------------------------
#  set.seed(123)
#  diaTrain <- biosign(dataMatrix[trainVi, ], sampleMetadata[trainVi, "type"],
#  bootI = 5)
#  set.seed(NULL)

## ----biosign_train, echo = FALSE, message = FALSE, warning = FALSE-------
set.seed(123)
diaTrain <- biosign(dataMatrix[trainVi, ], sampleMetadata[trainVi, "type"],
bootI = 5, plotL = FALSE)
set.seed(NULL)
plot(diaTrain, file.pdfC = "figures/diatrain.png")

## ----predict-------------------------------------------------------------
diaFitDF <- predict(diaTrain)

## ----confusion-----------------------------------------------------------
lapply(diaFitDF, function(predFc) table(actual = sampleMetadata[trainVi,
"type"], predicted = predFc))

## ----accuracy------------------------------------------------------------
sapply(diaFitDF, function(predFc) {
conf <- table(sampleMetadata[trainVi, "type"], predFc)
conf <- sweep(conf, 1, rowSums(conf), "/")
round(mean(diag(conf)), 3)
})

## ----getAccuracy---------------------------------------------------------
round(getAccuracyMN(diaTrain)["S", ], 3)

## ----performance---------------------------------------------------------
diaTestDF <- predict(diaTrain, newdata = dataMatrix[testVi, ])
sapply(diaTestDF, function(predFc) {
conf <- table(sampleMetadata[testVi, "type"], predFc)
conf <- sweep(conf, 1, rowSums(conf), "/")
round(mean(diag(conf)), 3)
})

## ----expressionset_code, eval = FALSE------------------------------------
#  library(Biobase)
#  diaSet <- ExpressionSet(assayData = t(dataMatrix),
#  phenoData = new("AnnotatedDataFrame", data = sampleMetadata))
#  set.seed(123)
#  biosign(diaSet, "type", bootI = 5)
#  set.seed(NULL)

## ----expressionset_figure, echo = FALSE, message = FALSE, warning = FALSE----
library(Biobase)
diaSet <- ExpressionSet(assayData = t(dataMatrix),
phenoData = new("AnnotatedDataFrame", data = sampleMetadata))
set.seed(123)
diaSet.bsg <- biosign(diaSet, "type", bootI = 5, plotL = FALSE)
set.seed(NULL)

## ----detach--------------------------------------------------------------
detach(diaplasma)

## ----sacurine_false, eval = FALSE----------------------------------------
#  data(sacurine)
#  set.seed(123) ##
#  sacSign <- biosign(sacurine[["dataMatrix"]],
#  sacurine[["sampleMetadata"]][, "gender"],
#  methodVc = "plsda")
#  set.seed(NULL)

## ----sacurine, echo = FALSE----------------------------------------------
data(sacurine)
set.seed(123) ##
sacSign <- biosign(sacurine[["dataMatrix"]],
sacurine[["sampleMetadata"]][, "gender"],
methodVc = "plsda", plotL = FALSE)
set.seed(NULL)
plot(sacSign, file.pdfC = "figures/sacsign.png")

## ----biomark, warning = FALSE, message = FALSE---------------------------
library(BioMark)
data(SpikePos)
group1Vi <- which(SpikePos[["classes"]] %in% c("control", "group1"))
appleMN <- SpikePos[["data"]][group1Vi, ]
spikeFc <- factor(SpikePos[["classes"]][group1Vi])
annotDF <- SpikePos[["annotation"]]
rownames(annotDF) <- colnames(appleMN)

## ----biomark_pca_false, eval = FALSE-------------------------------------
#  biomark.pca <- opls(appleMN, plotL = FALSE)
#  plot(biomark.pca, parAsColFcVn = spikeFc)

## ----biomark_pca, echo = FALSE-------------------------------------------
biomark.pca <- opls(appleMN, plotL = FALSE)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("overview", "outlier", "x-score", "x-loading"))
plot(biomark.pca, typeVc = typeC, parAsColFcVn = spikeFc, parDevNewL = FALSE)

## ----biomark_pls_false, eval = FALSE-------------------------------------
#  biomark.pls <- opls(appleMN, spikeFc)

## ----biomark_pls, echo = FALSE-------------------------------------------
biomark.pls <- opls(appleMN, spikeFc, plotL = FALSE)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("overview", "outlier", "x-score", "x-loading"))
plot(biomark.pls, typeVc = typeC, parDevNewL = FALSE)

## ----apple_biosign, eval = FALSE-----------------------------------------
#  set.seed(123)
#  appleSign <- biosign(appleMN, spikeFc)
#  set.seed(NULL)

## ----annotation, eval = FALSE--------------------------------------------
#  annotDF <- SpikePos[["annotation"]]
#  rownames(annotDF) <- colnames(appleMN)
#  annotDF[getSignatureLs(appleSign)[["complete"]], c("adduct", "found.in.standards")]

## ----golub_false, eval = FALSE-------------------------------------------
#  library(golubEsets)
#  data(Golub_Merge)
#  golubMN <- t(exprs(Golub_Merge))
#  leukemiaFc <- pData(Golub_Merge)[["ALL.AML"]]
#  table(leukemiaFc)
#  varSubVi <- 1501:2000
#  set.seed(123)
#  golubSign <- biosign(golubMN[, varSubVi], leukemiaFc, methodVc = "svm")
#  set.seed(NULL)

## ----golub, echo = FALSE, warning = FALSE, message = FALSE---------------
library(golubEsets)
data(Golub_Merge)
golubMN <- t(exprs(Golub_Merge))
leukemiaFc <- pData(Golub_Merge)[["ALL.AML"]]
table(leukemiaFc)
varSubVi <- 1501:2000
set.seed(123)
golubSign <- biosign(golubMN[, varSubVi], leukemiaFc, methodVc = "svm",
plotL = FALSE)
set.seed(NULL)
plot(golubSign, file.pdfC = "figures/golsign.png")

## ----hu6800, warning = FALSE, message = FALSE----------------------------
library(hu6800.db)
sapply(getSignatureLs(golubSign)[["complete"]],
       function(probeC)
       get(probeC, env = hu6800GENENAME))

## ----empty, echo = FALSE-------------------------------------------------
rm(list = ls())

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

