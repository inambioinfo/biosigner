test_biosign_plsda <- function(){

    data(diaplasma)

    varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))

    set.seed(123)

    biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                         factor(diaplasma[["sampleMetadata"]][, "type"]),
                         methodVc = "plsda",
                         bootI = 5)

    set.seed(NULL)

    plot(biosignLs)
    plot(biosignLs, tierMaxC = "A")
    plot(biosignLs, typeC = "boxplot")
    plot(biosignLs, tierMaxC = "A", typeC = "boxplot")

    checkEquals(biosignLs@tierMC["m427.215t07.9", "plsda"],
                "A")
    checkEqualsNumeric(biosignLs@accuracyMN["S", "plsda"],
                       0.7365702, tolerance = 1e-7)
    checkEqualsNumeric(getSummaryDF(biosignLs@modelLs[["plsda"]])[, "Q2(cum)"],
                       0.271, tolerance = 1e-6)

}

test_biosign_randomforest <- function(){

    data(diaplasma)

    varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))

    set.seed(123)

    biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                         factor(diaplasma[["sampleMetadata"]][, "type"]),
                         methodVc = "randomforest",
                         bootI = 5,
                         plotL = FALSE)

    set.seed(NULL)

    checkEquals(biosignLs@tierMC["m427.215t07.9", "randomforest"],
                "S")
    if (.Platform$OS.type != "windows") {
        checkEqualsNumeric(biosignLs@accuracyMN["AS", "randomforest"],
                           0.7999078, tolerance = 1e-7)
        checkEqualsNumeric(biosignLs@modelLs[["randomforest"]][["votes"]][2],
                           0.06470588, tolerance = 1e-6)
    }

}

test_biosign_svm <- function(){

    data(diaplasma)

    varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))

    set.seed(123)

    biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                         factor(diaplasma[["sampleMetadata"]][, "type"]),
                         methodVc = "svm",
                         bootI = 5)

    set.seed(NULL)

    checkEquals(biosignLs@tierMC["m123.998t01.0", "svm"],
                "A")
    checkEqualsNumeric(biosignLs@accuracyMN["AS", "svm"],
                       0.8313882, tolerance = 1e-7)
    checkEqualsNumeric(biosignLs@AS[["modelLs"]][["svm"]][["rho"]],
                       -0.7013267, tolerance = 1e-6)

}

test_biosign_predict <- function() {

    data(diaplasma)

    varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))

    samTotI <- nrow(diaplasma[["dataMatrix"]])
    trainVi <- 1:floor(samTotI/2)

    set.seed(123)

    biosignLs <- biosign(diaplasma[["dataMatrix"]][trainVi, varSelVi],
                         factor(diaplasma[["sampleMetadata"]][trainVi, "type"]),
                         bootI = 1,
                         plotL = FALSE)

    set.seed(NULL)

    predDF <- predict(biosignLs,
                      diaplasma[["dataMatrix"]][setdiff(1:samTotI, trainVi), varSelVi])
    predDIA043Vc <- as.character(unlist(predDF["DIA043", ]))
    if (.Platform$OS.type != "windows") {
        checkEquals(predDIA043Vc,
                    c("T2", "T1"))
    }

}


