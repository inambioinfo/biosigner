test_biosign_plsda <- function(){

    winL <- FALSE ## unit tests silenced on windows platforms because of errors on the moscato2 bioc platform running on windows 8

    data(diaplasma)

    varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))

    set.seed(123)

    biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                         diaplasma[["sampleMetadata"]][, "type"],
                         methodVc = "plsda",
                         bootI = 5)

    set.seed(NULL)

    plot(biosignLs)
    plot(biosignLs, tierMaxC = "A")
    plot(biosignLs, typeC = "boxplot")
    plot(biosignLs, tierMaxC = "A", typeC = "boxplot")


    if(winL || .Platform$OS.type != "windows") {

        checkEquals(biosignLs@tierMC["m427.215t07.9", "plsda"],
                    "A")
        checkEqualsNumeric(biosignLs@accuracyMN["S", "plsda"],
                           0.7365702, tolerance = 1e-7)

        library(ropls)

        checkEqualsNumeric(getSummaryDF(biosignLs@modelLs[["plsda"]])[, "Q2(cum)"],
                           0.271, tolerance = 1e-6)

    }

}

test_biosign_randomforest <- function(){

    winL <- FALSE ## unit tests silenced on windows platforms because of errors on the moscato2 bioc platform running on windows 8

    data(diaplasma)

    varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))

    set.seed(123)

    biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                         diaplasma[["sampleMetadata"]][, "type"],
                         methodVc = "randomforest",
                         bootI = 5,
                         plotL = FALSE)

    set.seed(NULL)


    if(winL || .Platform$OS.type != "windows") {

        checkEquals(biosignLs@tierMC["m427.215t07.9", "randomforest"],
                    "S")

        checkEqualsNumeric(biosignLs@accuracyMN["AS", "randomforest"],
                           0.7999078, tolerance = 1e-7)
        checkEqualsNumeric(biosignLs@modelLs[["randomforest"]][["votes"]][2],
                           0.06470588, tolerance = 1e-6)

    }

}

test_biosign_svm <- function(){

    winL <- FALSE ## unit tests silenced on windows platforms because of errors on the moscato2 bioc platform running on windows 8

    data(diaplasma)

    varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))

    set.seed(123)

    biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                         diaplasma[["sampleMetadata"]][, "type"],
                         methodVc = "svm",
                         bootI = 5)

    set.seed(NULL)


    if(winL || .Platform$OS.type != "windows") {

        checkEquals(biosignLs@tierMC["m123.998t01.0", "svm"],
                    "A")
        checkEqualsNumeric(biosignLs@accuracyMN["AS", "svm"],
                           0.8313882, tolerance = 1e-7)
        checkEqualsNumeric(biosignLs@AS[["modelLs"]][["svm"]][["rho"]],
                           -0.7013267, tolerance = 1e-6)

    }

}

test_biosign_predict <- function() {

    winL <- FALSE ## unit tests silenced on windows platforms because of errors on the moscato2 bioc platform running on windows 8

    data(diaplasma)

    varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))

    samTotI <- nrow(diaplasma[["dataMatrix"]])
    trainVi <- 1:floor(samTotI/2)

    set.seed(123)

    biosignLs <- biosign(diaplasma[["dataMatrix"]][trainVi, varSelVi],
                         diaplasma[["sampleMetadata"]][trainVi, "type"],
                         bootI = 1,
                         plotL = FALSE)

    set.seed(NULL)

    predDF <- predict(biosignLs,
                      diaplasma[["dataMatrix"]][setdiff(1:samTotI, trainVi), varSelVi])
    predDIA043Vc <- as.character(unlist(predDF["DIA043", ]))


    if(winL || .Platform$OS.type != "windows") {

        checkEquals(predDIA043Vc,
                    c("T2", "T1"))

    }

}

test_biosign_diaplasma <- function() {

    winL <- FALSE ## unit tests silenced on windows platforms because of errors on the moscato2 bioc platform running on windows 8

    data(diaplasma)

    varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))

    set.seed(123)

    biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                         diaplasma[["sampleMetadata"]][, "type"],
                         bootI = 5)

    set.seed(NULL)

    sigLs <- list(plsda = "m189.040t01.2",
                  randomforest = c("m427.215t07.9", "m189.040t01.2", "m995.613t10.2", "m455.221t08.1"),
                  svm = "m427.215t07.9",
                  complete = c("m427.215t07.9", "m189.040t01.2", "m995.613t10.2", "m455.221t08.1"))


    if(winL || .Platform$OS.type != "windows") {

        checkIdentical(getSignatureLs(biosignLs), sigLs)

        accMC <- matrix(c("0.7234427", "0.788674", "0.7384875", "0.7143228", "0.7566736", "0.8384932", "0.709905", "0.8271753", "0.6960716"),
                        nrow = 3,
                        ncol = 3,
                        dimnames = list(c("Full", "AS", "S"),
                            c("plsda", "randomforest", "svm")))

        biosignMC <- round(getAccuracyMN(biosignLs), 7)
        mode(biosignMC) <- "character"



        checkIdentical(biosignMC, accMC)

    }

}

test_biosign_sacurine <- function() {

    winL <- FALSE ## unit tests silenced on windows platforms because of errors on the moscato2 bioc platform running on windows 8

    library(ropls)

    data(sacurine)

    set.seed(123)

    biosignLs <- biosign(sacurine[["dataMatrix"]],
                         sacurine[["sampleMetadata"]][, "gender"],
                         bootI = 5)

    set.seed(NULL)

    sigLs <- list(plsda = c("Oxoglutaric acid", "Testosterone glucuronide", "p-Anisic acid", "Pantothenic acid", "Acetylphenylalanine", "Malic acid", "alpha-N-Phenylacetyl-glutamine", "Citric acid", "Gluconic acid and/or isomers", "Glucuronic acid and/or isomers", "Hippuric acid", "Phe-Tyr-Asp (and isomers)", "Threonic acid/Erythronic acid"),
                  randomforest = c("Oxoglutaric acid", "Testosterone glucuronide"),
                  svm = c("Oxoglutaric acid", "Testosterone glucuronide", "p-Anisic acid"),
                  complete = c("Oxoglutaric acid", "Testosterone glucuronide", "p-Anisic acid", "Pantothenic acid", "Acetylphenylalanine", "Malic acid", "alpha-N-Phenylacetyl-glutamine", "Citric acid", "Gluconic acid and/or isomers", "Glucuronic acid and/or isomers", "Hippuric acid", "Phe-Tyr-Asp (and isomers)", "Threonic acid/Erythronic acid"))


    if(winL || .Platform$OS.type != "windows") {

        checkIdentical(getSignatureLs(biosignLs), sigLs)

        accMC <- matrix(c("0.8958831", "0.8742228", "0.8745964", "0.847874", "0.8923372", "0.8482152", "0.8967789", "0.8667051", "0.8635216"),
                        nrow = 3,
                        ncol = 3,
                        dimnames = list(c("Full", "AS", "S"),
                            c("plsda", "randomforest", "svm")))

        biosignMC <- round(getAccuracyMN(biosignLs), 7)
        mode(biosignMC) <- "character"

        checkIdentical(biosignMC, accMC)

    }

}


