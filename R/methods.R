#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "ExpressionSet"),
          function(x, y, ...) {

              datMN <- t(exprs(x))

              samDF <- pData(x)

              if(!(y %in% colnames(samDF))) {
                  stop("'y' must be the name of a column of the phenoData slot of the 'ExpressionSet' object", call. = FALSE)
              } else {
                  rspFcVcn <- samDF[, y]
                  bsg <- biosign(datMN, rspFcVcn, ...)
              }

              bsg

          })

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "data.frame"),
          function(x, ...) {
              if(!all(sapply(x, data.class) == "numeric")) {
                  stop("'x' data frame must contain columns of 'numeric' vectors only", call. = FALSE)
              } else
                  x <- as.matrix(x)
              bsg <- biosign(x, ...)
              bsg
          })


#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "matrix"),
function(x,
         y,
         methodVc = c("all", "plsda", "randomforest", "svm")[1],
         bootI = 50,
         pvalN = 0.05,

         permI = 1,
         fixRankL = FALSE,

         printL = TRUE,
         plotL = TRUE,

         .sinkC = NULL,
         ...) {

    if(!is.null(.sinkC)) ##  Diversion of messages is required for the integration into Galaxy
        sink(.sinkC, append = TRUE)

    if(mode(x) != "numeric")
        stop("'x' matrix must be of 'numeric' mode", call. = FALSE)
    xMN <- x

    if(any(apply(xMN, 2, function(colVn) all(is.na(colVn)))))
        stop("'x' contains columns with 'NA' only", call. = FALSE)

    if(is.null(colnames(xMN))) {
        colnames(xMN) <- paste0("V", 1:ncol(xMN))
        warning(paste0("Missing column names in 'x' (i.e., variable names): Names are being set to '",
                       head(colnames(xMN), 1),
                       "' ... '",
                       tail(colnames(xMN), 1),
                       "'"),
                call. = FALSE)
    }

    if(is.character(y)) {
        yFc <- factor(y)
        warning("'y' character vector converted to a factor with levels '",
                paste(levels(y), collapse = "', '"),
                "'",
                call. = FALSE)
    }

    if(!is.factor(y)) {
        stop("'y' must be a factor", call. = FALSE)
    } else if(length(y) != nrow(xMN)) {
        stop("'y' factor length must be equal to the number of rows of 'x'", call. = FALSE)
    } else if(length(levels(y)) != 2) {
        stop("'y' must have two levels", call. = FALSE)
    } else if(any(is.na(y))) {
        stop("'y' must not contain missing ('NA') values")
    } else
        yFc <- y

    methodAllVc <- c("plsda", "randomforest", "svm")
    if(length(methodVc) == 1 && methodVc == "all") {
        methodVc <- methodAllVc
    } else if(!all(methodVc %in% methodAllVc))
        stop("The following method(s) is/are not available: '",
             paste(methodVc[!(methodVc %in% methodAllVc)], collapse = "', '"),
             "'",
             call. = FALSE)

    if(pvalN < 0 || pvalN > 1)
        stop("The p-value threshold 'pvalN' (",
             pvalN,
             ") must be between 0 and 1",
             call. = FALSE)


    ## signatures

    datasetLs <- list(dataMatrix = xMN,
                      sampleMetadata = data.frame(row.names = rownames(xMN),
                          respC = yFc),
                      variableMetadata = data.frame(row.names = colnames(xMN),
                          name. = colnames(xMN)))
    tierMN <- matrix(0, nrow = ncol(xMN), ncol = length(methodVc))
    rownames(tierMN) <- colnames(xMN)
    colnames(tierMN) <- methodVc

    accuVn <- numeric(length(methodVc))
    stopIterVl <- logical(length(methodVc))
    names(stopIterVl) <- names(accuVn) <- methodVc

    for(methodC in methodVc) {

        tierMetLs <- getTierF(datasetLs = datasetLs,
                              methodC = methodC,
                              bootI = bootI,
                              permI = permI,
                              pvalN = pvalN,
                              fixRankL = fixRankL)

        tierMN[, methodC] <- tierMetLs[["tierVn"]]
        accuVn[methodC] <- tierMetLs[["accuracyN"]]
        stopIterVl[methodC] <- tierMetLs[["stopIterL"]]

    }

    ## Formatting, ordering, and displaying results
    ##---------------------------------------------

    xSubMNAS <- xSubMN <- signatureLsAS <- signatureLs <- modelLsAS <- modelLs <- tierMC <- NULL

    ## Accuracy of full, S+A and S models

    accuracyMN <- matrix(NA, nrow = 3, ncol = length(methodVc),
                         dimnames = list(c("Full", "AS", "S"), methodVc))

    modelLs <- vector(mode = "list", length = length(methodVc))
    names(modelLs) <- methodVc
    signatureLs <- modelLsAS <- modelLs ## void models are empty lists
    for(sgnI in 1:length(signatureLs))
        signatureLs[[sgnI]] <- character() ## void signatures are character(0)
    signatureLsAS <- signatureLs

    ## translating tiers (building tierMC)

    if(max(tierMN) > 0) {
        for(methodC in methodVc) {
            tierVn <- tierMN[, methodC]
            translate <- rep("E", length(unique(tierVn)))
            if(max(tierVn) > 0){
                ## if(stop.rec[[methodC]]){
                if(stopIterVl[methodC]){
                    translate[max(tierVn)+1] <- "S"
                    if(max(tierVn) > 1)
                        translate[max(tierVn)] <- "A"
                    if(max(tierVn) > 2)
                        translate[max(tierVn)-1] <- "B"
                    if(max(tierVn) > 3)
                        translate[max(tierVn)-2] <- "C"
                    if(max(tierVn) > 4)
                        translate[max(tierVn)-3] <- "D"
                }
                else{
                    translate[max(tierVn)+1] <- "A"
                    if(max(tierVn) > 1)
                        translate[max(tierVn)] <- "B"
                    if(max(tierVn) > 2)
                        translate[max(tierVn)-1] <- "C"
                    if(max(tierVn) > 3)
                        translate[max(tierVn)-2] <- "D"
                }
            }
            tierVn <- translate[tierVn + 1]
            tierMC <- cbind(tierMC, tierVn)
        }
        dimnames(tierMC) <- dimnames(tierMN)

        ## ordering tiers (ordering tierMC)

        tierFullVc <- c("S", LETTERS[1:5])
        tierFullVn <- length(tierFullVc):1
        names(tierFullVn) <- tierFullVc
        tierOrdMC <- NULL
        for(tierC in tierFullVc) {
            if(any(tierC %in% c(tierMC))) {

                rowSelVl <- rowSums(tierMC == tierC) > 0
                tierMaxMC <- tierMC[rowSelVl, , drop = FALSE]
                tierMaxMN <- tierMaxMC
                for(j in 1:ncol(tierMaxMC))
                    tierMaxMN[, j] <- tierFullVn[tierMaxMC[, j]]
                mode(tierMaxMN) <- "numeric"
                tierMaxMC <- tierMaxMC[order(rowSums(tierMaxMN), decreasing = TRUE), , drop = FALSE]
                tierOrdMC <- rbind(tierOrdMC, tierMaxMC)
                tierMC <- tierMC[!rowSelVl, , drop = FALSE]
            }
        }
        tierMC <- tierOrdMC

        rm(tierMN)

        ## Accuracy of full, S+A and S models

        for(methodC in methodVc) {
            ## accuracyMN["Full", methodC] <- fsiResLs[[methodC]]$accuracyN
            accuracyMN["Full", methodC] <- accuVn[methodC]

            ## 'S+A' outputs
            datasetMethodLs <- datasetLs

            signatureLsAS[[methodC]] <- rownames(tierMC)[which(tierMC[, methodC] %in% c("S","A"))]

            if(length(signatureLsAS[[methodC]])) {
                datasetMethodLs$variableMetadata <- datasetMethodLs$variableMetadata[signatureLsAS[[methodC]], , drop = FALSE]
                datasetMethodLs$dataMatrix  <- datasetMethodLs$dataMatrix[ , signatureLsAS[[methodC]], drop = FALSE]
                my.fsi <- getBootSignificanceF(datasetMethodLs,
                                           "respC",
                                           methodC,
                                           bootI = bootI,
                                           permI = 0,
                                           fullModelL = TRUE)
                accuracyMN["AS", methodC] <- my.fsi$accuracyN
                modelLsAS[[methodC]] <- my.fsi$model
            }

            ## 'S' outputs
            datasetMethodLs <- datasetLs

            signatureLs[[methodC]] <- rownames(tierMC)[which(tierMC[, methodC] == "S")]

            if(length(signatureLs[[methodC]])) {
                datasetMethodLs$variableMetadata <- datasetMethodLs$variableMetadata[signatureLs[[methodC]], , drop = FALSE]
                datasetMethodLs$dataMatrix  <- datasetMethodLs$dataMatrix[ , signatureLs[[methodC]], drop = FALSE]
                my.fsi <- getBootSignificanceF(datasetMethodLs,
                                           "respC",
                                           methodC,
                                           bootI = bootI,
                                           permI = 0,
                                           fullModelL = TRUE)
                accuracyMN["S", methodC] <- my.fsi$accuracyN
                modelLs[[methodC]] <- my.fsi$model
            }
        }
        signatureLs[["complete"]] <- rownames(tierMC)[apply(tierMC, 1, function(rowVc) sum(rowVc == "S") > 0)]
        xSubMN <- xMN[, signatureLs[["complete"]], drop = FALSE] ## for boxplotting
        signatureLsAS[["complete"]] <- rownames(tierMC)[apply(tierMC, 1, function(rowVc) sum(rowVc %in% c("S", "A")) > 0)]
        xSubMNAS <- xMN[, signatureLsAS[["complete"]], drop = FALSE] ## for boxplotting

    } else {

        tierMC <- matrix("E", nrow = nrow(tierMN), ncol = ncol(tierMN),
                         dimnames = dimnames(tierMN))

        ## Accuracy of full model

        for(methodC in methodVc) {

            accuracyMN["Full", methodC] <- accuVn[methodC]

            xSubMN <- xMN[, signatureLs[["complete"]], drop = FALSE] ## for boxplotting

            xSubMNAS <- xMN[, signatureLsAS[["complete"]], drop = FALSE] ## for boxplotting

        }

    }

    AS <- list(modelLs = modelLsAS,
               signatureLs = signatureLsAS,
               xSubMN = xSubMNAS)

    bsg <- new("biosign")
    bsg@methodVc <- methodVc
    bsg@accuracyMN <- accuracyMN
    bsg@tierMC <- tierMC
    bsg@yFc <- yFc
    bsg@modelLs <- modelLs
    bsg@signatureLs <- signatureLs
    bsg@xSubMN <- xSubMN
    bsg@AS <- AS

    ## Printing
    ##---------

    if(printL) {
        show(bsg)
        warnings()
    }

    ## Plotting
    ##---------

    ## if(!is.null(bsg@accuracyMN) && plotL)
    if(!all(is.na(bsg@accuracyMN["S", ])) && plotL)
        plot(bsg)

    ## Closing connection
    ##-------------------

    if(!is.null(.sinkC)) ## Used in the Galaxy module
        sink()

    ## Returning
    ##----------

    return(invisible(bsg))


})


#' Plot method for 'biosign' signature objects
#'
#' Displays classifier tiers or individual boxplots from selected features
#'
#' @aliases plot.biosign plot,biosign-method
#' @param x An S4 object of class \code{biosign}, created by the \code{biosign}
#' function.
#' @param y Currently not used.
#' @param tierMaxC Character: Maximum level of tiers to display: Either 'S' and
#' 'A', (for boxplot), or also 'B', 'C', 'D', and 'E' (for tiers) by decreasing
#' number of selections
#' @param typeC Character: Plot type; either 'tier' [default] displaying the
#' comparison of signatures up to the selected 'tierMaxC' or 'boxplot' showing
#' the individual boxplots of the features selected by all the classifiers
#' @param file.pdfC Character: Figure filename ending with '.pdf'; default is
#' NULL (no saving; displaying instead)
#' @param .sinkC Character: Name of the file for R output diversion [default =
#' NULL: no diversion]; Diversion of messages is required for the integration
#' into Galaxy
#' @param ... Currently not used.
#' @return A plot is created on the current graphics device.
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @examples
#'
#' ## loading the diaplasma dataset
#'
#' data(diaplasma)
#' attach(diaplasma)
#'
#' ## restricting to a smaller dataset for this example
#'
#' featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
#' dataMatrix <- dataMatrix[, featureSelVl]
#' variableMetadata <- variableMetadata[featureSelVl, ]
#'
#' ## signature selection for all 3 classifiers
#' ## a bootI = 5 number of bootstraps is used for this example
#' ## we recommend to keep the default bootI = 50 value for your analyzes
#'
#' set.seed(123)
#' diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)
#'
#' ## individual boxplot of the selected signatures
#'
#' plot(diaSign, typeC = "boxplot")
#'
#' detach(diaplasma)
#'
#' @rdname plot
#' @export
setMethod("plot", signature(x = "biosign"),
          function(x,
                   y,
                   tierMaxC = "S",
                   typeC = c("tier", "boxplot")[1],
                   file.pdfC = NULL,
                   .sinkC = NULL,
                   ...) {

    if(!is.null(.sinkC)) ##  Diversion of messages is required for the integration into Galaxy
        sink(.sinkC, append = TRUE)

    tierFullVc <- c("S", LETTERS[1:5])

    if(length(tierMaxC) != 1 || any(!(tierMaxC %in% tierFullVc))) {
        message("'tierMaxC' argument must be either '", paste(tierFullVc, collapse = "', '"), "' for the 'tier' plot")
        if(!is.null(.sinkC)) ## Used in the Galaxy module
            sink()
        return(invisible(NULL))
    } else if(typeC == "boxplot" && any(!(tierMaxC %in% c("S", "A")))) {
        message("'tierMaxC' argument must be either 'S' or 'A' for the 'boxplot'")
        if(!is.null(.sinkC)) ## Used in the Galaxy module
            sink()
        return(invisible(NULL))
    } else
        tierVc <- tierFullVc[1:which(tierFullVc == tierMaxC)]

    switch(typeC,

           tier = {

               if(sum(x@tierMC %in% tierVc) == 0) {

                   message("No signature up to tier '", tierMaxC, "' to be plotted")
                   if(!is.null(.sinkC)) ## Used in the Galaxy module
                       sink()
                   return(invisible(NULL))

               }

           },
           boxplot = {

               switch(tierMaxC,
                      S = {
                          sgnAllVc = x@signatureLs[["complete"]]
                          xSubMN <- x@xSubMN
                      },
                      A = {
                          sgnAllVc = x@AS[["signatureLs"]][["complete"]]
                          xSubMN <- x@AS[["xSubMN"]]
                        })

               if(length(sgnAllVc) == 0) {
                   message("No signature up to tier '", tail(tierVc, 1), "' to be plotted")
                   if(!is.null(.sinkC)) ## Used in the Galaxy module
                       sink()
                   return(invisible(NULL))
               }

           })


    if(is.null(file.pdfC)) {
        dev.new()
    } else {
        if(grepl("pdf$", file.pdfC)) {
            pdf(file.pdfC)
        } else if(grepl("png$", file.pdfC)) {
            png(file.pdfC)
        } else
            stop("'file.pdfC' argument must have a '.pdf' extension",
                 call. = FALSE)
    }

    switch(typeC,

           tier = {

               if(length(setdiff(tierVc, c(x@tierMC)))) {
                   tierNotFoundVc <- tierFullVc[tierFullVc %in% setdiff(tierVc, c(x@tierMC))]
                   warning("tierMC does not contain the following values: '", paste(tierNotFoundVc, collapse = "', '"), "'", call. = FALSE)
               }
               inputMC <- x@tierMC[apply(x@tierMC, 1, function(rowVc) sum(rowVc %in% tierVc) > 0), ,
                                        drop = FALSE]
               tierFullVi <- 1:length(tierFullVc)
               names(tierFullVi) <- tierFullVc
               inputMN <- inputMC
               for(j in 1:ncol(inputMN))
                   inputMN[, j] <- tierFullVi[inputMC[, j]]
               mode(inputMN) <- "numeric"


               imageMN <- inputMN

               colnames(imageMN) <- gsub("plsda", "PLS-DA",
                                         gsub("randomforest", "Random Forest",
                                              gsub("svm", "SVM", colnames(imageMN))))

               if(length(rownames(imageMN)) == 0) {

                   rownameInVc <- rep("", times = nrow(imageMN))

               } else {

                   rownameInVc <- rownames(imageMN)

                   warnOpN <- getOption("warn")
                   options(warn = -1)
                   rownameInCharVsNumL <- NA %in% as.numeric(rownameInVc)
                   options(warn = warnOpN)

                   if(rownameInCharVsNumL) {

                       rownameInDuplicateVc <- duplicated(rownameInVc)

                       rownameInVc <- sapply(1:length(rownameInVc), function(i) ifelse(!rownameInDuplicateVc[i],
                                                                                       yes = rownameInVc[i],
                                                                                       no = ""))

                   }
               }

               colnameInVc <- colnames(imageMN)


               imageMN <- t(imageMN[nrow(imageMN):1, , drop = FALSE])

               par(bg = "white",
                   font = 2,
                   lwd  = 2)

               layout(matrix(c(2, 1),
                             nrow = 1,
                             ncol = 2,
                             byrow = TRUE),
                      widths = c(8, 2))

               tierFullColVc <- c("#1A9850", "#91CF60", "#D9EF8B", "#FEE08B", "#FC8D59", "#D73027")
               names(tierFullColVc) <- tierFullVc

               ## draw color scale

               par(mar = c(1.1, 0.6, 7.6, 4.1))

               plot(x = 0,
                    y = 0,
                    font.axis = 2,
                    font.lab = 2,
                    type = "n",
                    xlim = c(0, 1),
                    ylim = c(0, 6),
                    xlab = "",
                    ylab = "",
                    xaxs = "i",
                    yaxs = "i",
                    xaxt = "n",
                    yaxt = "n")

               rect(xleft = 0,
                    ybottom = 0:5,
                    xright = 1,
                    ytop = 1:6,
                    col = rev(tierFullColVc),
                    border = NA)

               axis(at = tierFullVi - 0.5,
                    font = 2,
                    font.axis = 2,
                    labels = rev(tierFullVc),
                    las = 1,
                    lwd = 2,
                    lwd.ticks = 2,
                    side = 4)

               arrows(par("usr")[2],
                      par("usr")[4],
                      par("usr")[2],
                      par("usr")[3],
                      code = 0,
                      lwd = 2,
                      xpd = TRUE)

               ## draw image

               par(mar = c(1.1,
                       17,
                       7.6,
                       0.3))

               image(x = 1:nrow(imageMN),
                     y = 1:ncol(imageMN),
                     z = imageMN,
                     col = tierFullColVc[min(inputMN):max(inputMN)],
                     font.axis = 2,
                     font.lab = 2,
                     xaxt = "n",
                     yaxt = "n",
                     xlab = "",
                     ylab = "")

               xLabelAtVn <- 1:ncol(imageMN)
               yLabelAtVn <- 1:nrow(imageMN)

               ## if(colnameInCharVsNumL) {

               xLabelVc <- colnameInVc[colnameInVc != ""]

               xLabelLowIndiceVn <- which(colnameInVc != "")

               xLabelSpanVn <- diff(c(xLabelLowIndiceVn, length(colnameInVc) + 1))

               xLabelAtVn <- xLabelLowIndiceVn - rep(1, times = length(xLabelLowIndiceVn)) + xLabelSpanVn / 2 + rep(0.5, times = length(xLabelLowIndiceVn))

               par(cex = 1)

               axis(side = 3,
                    at = xLabelAtVn,
                    font = 2,
                    labels = xLabelVc,
                    las = 2,
                    line = -0.5,
                    tick = FALSE)

               par(cex = 1)


               if(rownameInCharVsNumL) {

                   yLabelVc <- rownameInVc[rownameInVc != ""]

                   yLabelLowIndiceVn <- which(rownameInVc != "")

                   yLabelSpanVn <- diff(c(yLabelLowIndiceVn, length(rownameInVc) + 1))

                   yLabelAtVn <- ncol(imageMN) - rev(yLabelLowIndiceVn) + rep(1, times = length(yLabelLowIndiceVn)) - rev(yLabelSpanVn) / 2 + rep(0.5, times = length(yLabelLowIndiceVn))

                   par(cex = 1)

                   axis(side = 2,
                        at = yLabelAtVn,
                        font = 2,
                        hadj = 1,
                        labels = rev(yLabelVc),
                        las = 2,
                        line = -0.5,
                        tick = FALSE)

               } else {

                   rownameVn <- as.numeric(rownameInVc)

                   prettyVn <- pretty(rownameVn)

                   prettyVn <- prettyVn[min(rownameVn) <= prettyVn & prettyVn <= max(rownameVn)]

                   indiceVn <- numeric()

                   for(k in 1:length(prettyVn))
                       indiceVn[k] <- which(abs(rownameVn - prettyVn[k]) == min(abs(rownameVn - prettyVn[k])))[1]

                   axis(side = 2,
                        at = nrow(imageMN) - rev(indiceVn),
                        font = 2,
                        labels = as.character(prettyVn))

               }

               par(cex = 1)

               ## additional lines

               if(rownameInCharVsNumL)
                   abline(h = ncol(imageMN) - which(rownameInVc != "") + 1 + 0.5,
                          lwd = 1)

               ## if(colnameInCharVsNumL)
               abline(v = which(colnameInVc != "") - 0.5,
                      lwd = 1)


               ## border

               box(lwd = 2)

               ## arrows at the end of the axes

               arrows(par("usr")[1],
                      par("usr")[4],
                      par("usr")[1],
                      par("usr")[3],
                      length = 0.1,
                      lwd = 2,
                      xpd = TRUE)

               box(lwd=3)

           },
           boxplot = {

               palVc <- c(code100_plsda = "#CB181D",
                          ## code110_plsda_randomforest = "#F16913",
                          code010_randomforest = "#238B45",
                          ## code011_randomforest_svm = "#807DBA",
                          code001_svm = "#2171B5")
                          ## code101_plsda_svm = "          #BF812D",
                          ## code111_plsda_randomforest_svm = "black"
                          ## )

               palDF <- data.frame(row.names = substr(names(palVc), 5, 7),
                                   classifier = substr(names(palVc), 9, nchar(names(palVc))),
                                   color = palVc,
                                   stringsAsFactors = FALSE)

               dimN <- ceiling(sqrt(length(sgnAllVc)))
               layout(matrix(1:dimN^2, nrow = dimN, byrow = TRUE))
               par(font = 2, font.axis = 2, font.lab = 2, las = 1,
                   mar = c(2.1, 2.6, 2.6, 1.1),
                   oma = c(0, 0, 2.1, 0))
               for(bmkC in sgnAllVc) {
                   bmkTirVc <- x@tierMC[bmkC, ]
                   names(bmkTirVc) <- colnames(x@tierMC)
                   bmkTirSigVl <- bmkTirVc %in% tierVc
                   bmkModVc <- paste(names(bmkTirVc)[bmkTirSigVl], collapse = ", ")
                   bmkModVc <- gsub("plsda", "PLSDA", gsub("randomforest", 'RF', gsub("svm", "SVM", bmkModVc)))
                   bmkColC <- ifelse(sum(bmkTirSigVl) == 1,
                                     palDF[paste(as.numeric(bmkTirSigVl), collapse = ""), "color"],
                                     "black")
                   boxplot(xSubMN[, bmkC] ~ x@yFc,
                           border = bmkColC,
                           main = "")
                   mtext(ifelse(nchar(bmkC) > 23, paste0(substr(bmkC, 1, 20), "."), bmkC),
                         cex = 0.8,
                         line = 1.1)
                   mtext(bmkModVc, line = 0.1, cex = 0.6)
               }
               title(main = paste0("'", switch(tierMaxC, S = "S", A = "S+A"), "' signature"), line = 0.5, cex.main = 1.5, outer = TRUE)

           })


    if(!is.null(file.pdfC))
        dev.off()


    ## Closing connection
    ##-------------------

    if(!is.null(.sinkC)) ## Used in the Galaxy module
        sink()

          })


#' Predict method for 'biosign' signature objects
#'
#' This function predicts values based upon \code{biosign} classifiers trained
#' on the 'S' signature
#'
#'
#' @aliases predict.biosign predict,biosign-method
#' @param object An S4 object of class \code{biosign}, created by
#' \code{biosign} function.
#' @param newdata Either a data frame or a matrix, containing numeric columns
#' only, with column names identical to the 'x' used for model training with
#' 'biosign'.
#' @param tierMaxC Character: Maximum level of tiers to display: Either 'S'or
#' 'A'.
#' @param ... Currently not used.
#' @return Data frame with the predictions for each classifier as factor
#' columns.
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @examples
#'
#' ## loading the diaplasma dataset
#'
#' data(diaplasma)
#' attach(diaplasma)
#'
#' ## restricting to a smaller dataset for this example
#'
#' featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
#' dataMatrix <- dataMatrix[, featureSelVl]
#' variableMetadata <- variableMetadata[featureSelVl, ]
#'
#' ## training the classifiers
#' ## a bootI = 5 number of bootstraps is used for this example
#' ## we recommend to keep the default bootI = 50 value for your analyzes
#'
#' set.seed(123)
#' diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)
#'
#' ## fitted values (for the subsets restricted to the 'S' signatures)
#' sFitDF <- predict(diaSign)
#'
#' ## confusion tables
#' print(lapply(sFitDF, function(predFc) table(actual = sampleMetadata[,
#' "type"], predicted = predFc)))
#'
#' ## balanced accuracies
#' sapply(sFitDF, function(predFc) { conf <- table(sampleMetadata[,
#' "type"], predFc)
#' conf <- sweep(conf, 1, rowSums(conf), "/")
#' mean(diag(conf))
#' })
#' ## note that these values are slightly different from the accuracies
#' ## returned by biosign because the latter are computed by using the
#' ## resampling scheme selected by the bootI or crossvalI arguments
#' getAccuracyMN(diaSign)["S", ]
#'
#' detach(diaplasma)
#'
#' @rdname predict
#' @export
setMethod("predict", signature(object = "biosign"),
          function(object, newdata, tierMaxC = "S", ...)
          {

    if(any(!(tierMaxC %in% c("S", "A"))))
        stop("'tierMaxC' argument must be 'S' or 'A' for predictions", call. = FALSE)

    switch(tierMaxC,
           S = {
               signatureLs <- object@signatureLs
               modelLs <- object@modelLs
           },
           A = {
               signatureLs <- object@AS[["signatureLs"]]
               modelLs <- object@AS[["modelLs"]]
           })

    if(missing(newdata)) {

        fitLs <- lapply(modelLs,
                        function(model) {
                            if(is.null(model))
                                return(NULL)
                            else
                                return(switch(class(model),
                                              opls = fitted(model),
                                              randomForest = model[["predicted"]],
                                              svm = model[["fitted"]]))
                        })

        fitLs <- fitLs[!sapply(fitLs, is.null)]

        if(length(fitLs) == 0) {
            warning("Empty signatures for all classifiers up to tier '", tierMaxC, "'; fitted output is set to NULL", call. = FALSE)
            return(invisible(NULL))
        } else
            return(as.data.frame(fitLs))

    } else {

        if(is.data.frame(newdata)) {
            if(!all(sapply(newdata, data.class) == "numeric")) {
                stop("'newdata' data frame must contain numeric columns only", call. = FALSE)
            }
        } else if(is.matrix(newdata)) {
            if(mode(newdata) != "numeric") {
                stop("'newdata' matrix must be of 'numeric' mode", call. = FALSE)
            } else
                newdata <- as.data.frame(newdata)
        } else
            stop("'newdata' must be either a data.frame or a matrix", call. = FALSE)

        if(is.null(colnames(newdata))) {
            colnames(newdata) <- paste0("V", 1:ncol(newdata))
            warning(paste0("Missing column names in 'newdata' (i.e., variable names): Names are being set to '",
                           head(colnames(newdata), 1),
                           "' ... '",
                           tail(colnames(newdata), 1),
                           "'"),
                    call. = FALSE)
        }

        if(length(signatureLs[["complete"]]) == 0){

            warning("Signatures from all classifiers up to tier '", tierMaxC, "' are empty; prediction output is set to NULL", call. = FALSE)

            return(invisible(NULL))

        } else {

            predDF <- as.data.frame(lapply(object@methodVc,
                                           function(methodC) {
                                               signVc <- signatureLs[[methodC]]
                                               if(length(signVc)) {
                                                   signOutVc <- signVc[!(signVc %in% colnames(newdata))]
                                                   if(length(signOutVc)) {
                                                       stop("The following variables from the ", methodC, " ", switch(tierMaxC, S = "S", A = "S+A"), " signature are not found as column names in 'newdata': '", paste(signOutVc, collapse = "', '"), "'.", call. = FALSE)
                                                   } else
                                                       return(predict(modelLs[[methodC]],
                                                                      newdata = newdata[, signVc, drop = FALSE]))
                                               } else
                                                   return(rep(NA, nrow(newdata))) ## /!\ logical mode
                                           }))
            colnames(predDF) <- object@methodVc
            predDF <- predDF[, !apply(predDF, 2, function(predVn) all(is.na(predVn))), drop = FALSE]

            return(predDF)

        }

    }

          })


#' Show method for 'biosign' signature objects
#'
#' Prints the selected features and the accuracies of the classifiers.
#'
#' @aliases show.biosign show,biosign-method
#' @param object An S4 object of class \code{biosign}, created by the
#' \code{biosign} function.
#' @return Invisible.
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @examples
#'
#' ## loading the diaplasma dataset
#'
#' data(diaplasma)
#' attach(diaplasma)
#'
#' ## restricting to a smaller dataset for this example
#'
#' featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
#' dataMatrix <- dataMatrix[, featureSelVl]
#' variableMetadata <- variableMetadata[featureSelVl, ]
#'
#' ## signature selection for all 3 classifiers
#' ## a bootI = 5 number of bootstraps is used for this example
#' ## we recommend to keep the default bootI = 50 value for your analyzes
#'
#' set.seed(123)
#' diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)
#'
#' diaSign
#'
#' detach(diaplasma)
#'
#' @rdname show
#' @export
setMethod("show",
          "biosign",
          function(object)
          {

              tierMaxC <- "S"

              ## if(is.null(object@accuracyMN)) {
              if(all(is.na(object@accuracyMN["S", ]))) {

                  cat("No significant variable found for the selected classifier(s): '", paste(object@methodVc, collapse = "', '"), "'\n", sep = "")

              } else {

                  tierFullVc <- c("S", "A", "B", "C", "D", "E")

                  if(any(!(tierMaxC %in% tierFullVc))) {
                      stop("'tierMaxC' argument must be in '", paste(tierFullVc, collapse = "', '"), "'", call. = FALSE)
                  } else
                      tierVc <- tierFullVc[1:which(tierFullVc == tierMaxC)]

                  if(sum(object@tierMC %in% tierVc)) {
                      if(length(setdiff(tierVc, c(object@tierMC)))) {
                          tierNotFoundVc <- tierFullVc[tierFullVc %in% setdiff(tierVc, c(object@tierMC))]
                          warning("tierMC does not contain the following values: '", paste(tierNotFoundVc, collapse = "', '"), "'", call. = FALSE)
                      }

                      cat("Significant features from '", paste(tierVc, collapse = "', '"), "' groups:\n", sep = "")
                      print(object@tierMC[apply(object@tierMC, 1, function(rowVc) sum(rowVc %in% tierVc) > 0), ,
                                          drop = FALSE])

                      cat("Accuracy:\n", sep = "")
                      print(round(object@accuracyMN, 3))

                      invisible(NULL)

                  }

              }

          })

#' @rdname getAccuracyMN
#' @export
setMethod("getAccuracyMN", "biosign",
          function(object) {
              return(object@accuracyMN)
          })

#' @rdname getSignatureLs
#' @export
setMethod("getSignatureLs", "biosign",
          function(object, tierC = c("S", "AS")[1]) {
              if(tierC == "S") {
                  return(object@signatureLs)
              } else if(tierC == "AS") {
                  return(object@AS$signatureLs)
              } else
                  stop("'tierC' argument must be either 'S' or 'AS'")
          })




