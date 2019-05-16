
seqgroups <- function(alldata = NULL, labels = NULL, lag = 1, adjacent = TRUE,
                      onezero = NULL, tailed = 2, test = "homogeneity", output = "all") {

cat('\n\nLag Sequential Analysis of data that are in segments (e.g, for multiple dyads or groups\n\n')

cat("\nRequested 'tail' (1 or 2) for Significance Tests =", tailed, "\n")

  if (!adjacent) {
    adjacent <- 0
  } else if (adjacent && is.null(onezero)) {
    adjacent <- 1
  } else if (adjacent && !is.null(onezero)) {
    adjacent <- 2
  }

  alldata <- matrix(alldata, ncol = 1)


  # are all data values numeric?
  if ((all(sapply(alldata, is.numeric))) == TRUE) {

    # create a version of alldata without the numeric segment codes, in order to get the max code value
    dattemp <- alldata[alldata < 1000]
    ncodes <- max(dattemp)

    if (is.null(labels)) {
      labels <- 1:ncodes
      for (lupe in 1:ncodes) labels[lupe] <- paste("Code", lupe)
    }
  }


  # if any data values are characters, treat them all as strings & provide integer values for the analyses
  if ((any(sapply(alldata, is.character))) == TRUE) {

    # create a version of alldata without the "segment" values, in order to get the code labels
    dattemp <- alldata[nzchar(gsub("segment+", "", alldata))]
    labels <- unique(dattemp)

    newdat <- matrix(-9999, length(alldata), 1)
    for (lupe in 1:length(alldata)) {
      # replace code values that aren't segment codes with integers
      if (alldata[lupe, 1] != "segment") newdat[lupe, 1] <- which(labels == alldata[lupe, 1], arr.ind = F)
      if (alldata[lupe, 1] == "segment") newdat[lupe, 1] <- 999999
    }
    alldata <- newdat
    ncodes <- length(labels)
  }




  # pooling data for stationarity analyses

  if (test == "stationarity") {
    totdata <- 0
    for (x in 1:nrow(alldata)) {
      if (alldata[x] < 1000) {
        totdata <- rbind(totdata, alldata[x])
      }
    }
    totdata <- matrix(totdata[2:nrow(totdata)], ncol = 1)
    freqst <- matrix(0, ncodes, ncodes)
    for (y in 1:length(totdata)) {
      if (y + lag <= length(totdata)) {
        freqst[totdata[y], totdata[y + lag]] <- freqst[totdata[y], totdata[y +
          lag]] + 1
      }
    }
  }


#  out <- matrix(0, 1, (ncodes * ncodes + 1))
  each <- matrix(0, 1, (ncodes * ncodes))
  start <- 2
  lastone <- 0
  for (a in 1:nrow(alldata)) {
    freqs <- matrix(0, ncodes, ncodes)
    data <- 0
    for (b in start:nrow(alldata)) {
      if (alldata[b] < 1000) {
        data <- rbind(data, alldata[b])
        end <- b
      }
      if (alldata[b] >= 1000) {
        start <- b + 1
        break
      }
    }
    data <- matrix(data[2:length(data)], ncol = 1)
    lastone <- rbind(lastone, data[nrow(data)])

    freqs <- matrix(0, ncodes, ncodes)
    for (c in 1:nrow(data)) {
      if (c + lag <= nrow(data)) {
        freqs[data[c], data[c + lag]] <- freqs[data[c], data[c + lag]] + 1
      }
    }
    each <- rbind(each, matrix(t(freqs), 1, (ncodes * ncodes)))
    if (end == nrow(alldata)) {
      break
    }
  }

  lastone <- matrix(c(lastone[2:nrow(lastone), 1], lastone[nrow(lastone), 1]), ncol = 1)



  if (test == "homogeneity") {
    freqst <- (matrix(t(colSums(each)), ncodes, ncodes))
    freqst <- t(freqst)
  }
  rowtott <- rowSums(freqst)
  pp <- matrix(0, ncodes, ncodes)
  each <- rbind(each[2:nrow(each), ], matrix(t(freqst), 1, ncol(each)))

  # initializing
  lrx2sh <- 0
  for (dindex in 1:nrow(each)) {
    lrx2t <- matrix(0)
    freqs <- t(matrix(t(each[dindex, ]), ncodes, ncodes))
    rowtots <- matrix(rowSums(freqs))
    coltots <- matrix(colSums(freqs), ncol = ncodes)
    ntrans <- sum(rowtots)
    prows <- rowtots / ntrans
    pcols <- coltots / ntrans
    tprob <- matrix(-9999, ncodes, ncodes)
    et <- matrix(-9999, ncodes, ncodes)
    expfreq <- matrix(-9999, ncodes, ncodes)
    zadjres <- matrix(-9999, ncodes, ncodes)
    pzadjres <- matrix(-1, ncodes, ncodes)
    yulesq <- matrix(-9999, ncodes, ncodes)
    var <- matrix(-9999, ncodes, ncodes)
    min <- matrix(-9999, ncodes, ncodes)
    kappa <- matrix(-9999, ncodes, ncodes)
    zkappa <- matrix(-9999, ncodes, ncodes)
    pzkappa <- matrix(1, ncodes, ncodes)
    n <- ntrans + 1
    nr <- rowtots
    nr[lastone[dindex, 1]] <- nr[lastone[dindex, 1]] + 1

    # Warning message for specification error.
    if ((adjacent == 0) & (any(diag(freqs) == 0))) {
      cat("\n\nWarning:\nYou have indicated that consecutive or adjacent codes can never repeat")
      cat("\n(adjacent = 0), yet repeating codes have been found in the data. See the ")
      cat("\nmain diagonal of the frequency matrix. This specification error will result ")
      cat("\nin faulty computations for LRX2, z-values, and adjusted residuals.\n")
    }

    for (i in 1:ncodes) {
      for (j in 1:ncodes) {

        # Note: more refined computations for when adjacent codes cannot repeat appear below,
        # after the above 2 loops are complete.
        if ((adjacent == 0) && (ntrans - rowtots[i] > 0)) {
          pcols[j] <- coltots[j] / (ntrans - rowtots[i])
        }
        if ((adjacent == 0) && (ntrans - rowtots[j] > 0)) {
          expfreq[i, j] <- (rowtots[i] * coltots[j]) / (ntrans - rowtots[j])
        }
        if ((adjacent == 0) && (n - nr[i] > 0)) {
          et[i, j] <- (nr[i] * nr[j]) / (n - nr[i])
        }

        if (adjacent == 1) {
          et[i, j] <- (nr[i] * nr[j]) / n
          expfreq[i, j] <- (rowtots[i] * coltots[j]) / ntrans
        }

        # transitional probabilities.
        if (rowtots[i] > 0) {
          tprob[i, j] <- freqs[i, j] / rowtots[i]
        }

        # tablewise LRX2
        if ((freqs[i, j] > 0) && (expfreq[i, j] > 0)) {
          lrx2t <- lrx2t + 2 * (freqs[i, j] * log(freqs[i, j] / expfreq[i, j]))
        }

        # LRX2 for stationarity/homogeneity.
        if (dindex < nrow(each)) {
          if (rowtott[i] > 0) {
            pp[i, j] <- freqst[i, j] / rowtott[i]
          }
          if ((tprob[i, j] > 0) && (pp[i, j] > 0)) {
            lrx2sh <- lrx2sh + 2 * (freqs[i, j] * log(tprob[i, j] / pp[i, j]))
          }
        }

        # adjusted residuals (z values) & sig levels.
        if ((expfreq[i, j] * (1 - pcols[j]) * (1 - prows[i])) > 0) {
          zadjres[i, j] <- ((freqs[i, j] - expfreq[i, j])) / sqrt(expfreq[i, j] *
            (1 - pcols[j]) * (1 - prows[i]))
          pzadjres[i, j] <- (1 - pnorm(abs(zadjres[i, j]))) * tailed
        }

        # Yule's Q; see Bakeman & Gottman, 1997, p. 129.
        a <- freqs[i, j]
        b <- rowtots[i] - freqs[i, j]
        c <- coltots[j] - freqs[i, j]
        d <- ntrans - rowtots[i] - coltots[j] + freqs[i, j]
        if ((a * d + b * c) > 0) {
          yulesq[i, j] <- (a * d - b * c) / (a * d + b * c)
        }

        # kappas, z values & sig levels.
        var[i, j] <- (nr[i] * nr[j] * (n - nr[j]) * (n - nr[i])) / (n^2 * (n - 1))
        if (var[i, j] > 0) {
          zkappa[i, j] <- (freqs[i, j] - et[i, j]) / sqrt(var[i, j])
          if (nr[i] <= nr[j]) {
            min[i, j] <- nr[i]
          } else {
            min[i, j] <- nr[j]
          }
          if (min[i, j] - et[i, j] != 0) {
            kappa[i, j] <- (freqs[i, j] - et[i, j]) / (min[i, j] - et[i, j])
            if (kappa[i, j] < 0) {
              kappa[i, j] <- (freqs[i, j] - et[i, j]) / et[i, j]
            }
            pzkappa[i, j] <- (1 - pnorm(abs(zkappa[i, j]))) * tailed
          }
        }
      }
    }

    if ((adjacent == 0) || (adjacent == 2)) {

      # maximum likelihood estimation of the expected cell frequencies using
      # iterative proportional fitting (Wickens, 1989, pp. 107-112).

      rsumsf <- rowSums(freqs)
      csumsf <- colSums(freqs)

      if (is.null(onezero)) {
        onezero <- matrix(1, ncodes, ncodes)
        diag(onezero) <- 0
      }

      # The two previous commands create a matrix of ones and
      # zeros that is used in estimating the expected cell
      # frequencies. A "1" indicates that the expected frequency
      # for a given cell is to be estimated, whereas a "0"
      # indicates that the expected frequency for the cell
      # should NOT be estimated, typically because it is a
      # structural zero (codes that cannot follow one another).
      # By default, the matrix that is created by the above
      # commands has zeros along the main diagonal, and ones
      # everywhere else, which will be appropriate for most data
      # sets. However, if yor data happen to involve structural
      # zeros that occur in cells other than the cells along the
      # main diagonal, then you must create a ONEZERO matrix
      # with ones and zeros that is appropriate for your data
      # before running any of the commands below. Enter your
      # ONEZERO matrix now.

      expfreq <- onezero

      for (ipfloop in 1:100) {

        # adjusting by row.
        xr <- matrix(0, ncodes, 1)
        rsumse <- rowSums(expfreq)
        for (r in 1:ncodes) {
          if (rsumse[r] > 0) {
            xr[r] <- rsumsf[r] / rsumse[r]
          }
        }
        for (i in 1:ncodes) {
          for (j in 1:ncodes) {
            if (onezero[i, j] == 1) {
              expfreq[i, j] <- expfreq[i, j] * xr[i]
            }
          }
        }

        # adjusting by column.
        xc <- matrix(0, 1, ncodes)
        csumse <- colSums(expfreq)
        for (c in 1:ncodes) {
          if (csumse[c] > 0) {
            xc[c] <- csumsf[c] / csumse[c]
          }
        }
        for (i in 1:ncodes) {
          for (j in 1:ncodes) {
            if (onezero[i, j] == 1) {
              expfreq[i, j] <- expfreq[i, j] * xc[j]
            }
          }
        }
        rdiffs <- rsumsf - rowSums(expfreq)
        cdiffs <- csumsf - colSums(expfreq)
        if ((max(rdiffs) < 1e-04) & (max(cdiffs) < 1e-04)) {
          break
        }
      }
      cat("Maximum likelihood estimation of the expected cell frequencies using iterative proportional fitting ")
      if ((max(rdiffs) < 1e-04) & (max(cdiffs) < 1e-04)) {
        cat("converged after the following number of iterations:", ipfloop, "\n\n")
      } else {
        cat(
          "did NOT converge after the following number of iterations:", ipfloop,
          "\n\n"
        )
      }

      # tablewise LRX2
      lrx2t <- matrix(0)
      for (i in 1:ncodes) {
        for (j in 1:ncodes) {
          if ((freqs[i, j] > 0) & (expfreq[i, j] > 0)) {
            lrx2t <- lrx2t + 2 * freqs[i, j] * log(freqs[i, j] / expfreq[i, j])
          }
        }
      }

      # adjusted residuals for matrices with structural zeros (Christensen, 1997, p. 357).

      # constructing the design matrix. 345
      x <- matrix(1, nrow = ncodes^2, ncol = 1)
      y <- matrix(0, nrow = ncodes^2, ncol = ncodes - 1)
      z <- matrix(0, nrow = ncodes^2, ncol = ncodes - 1)
      for (i in 1:(ncodes - 1)) {
        for (j in 1:ncodes) {
          y[i * ncodes + j, i] <- 1
          z[(((j - 1) * ncodes) + (i + 1)), i] <- 1
        }
      }
      des1 <- cbind(x, y, z)

      # pruning values corresponding to cells with structural zeros.
      onezero2 <- matrix(t(onezero), ncodes^2, 1)
      dm1 <- matrix(t(expfreq), ncodes^2, 1)
      dm2 <- matrix(-9999, 1, 1)
      des2 <- matrix(-9999, 1, ncol(des1))
      for (p in 1:(ncodes^2)) {
        if (onezero2[p] == 1) {
          dm2 <- rbind(dm2, dm1[p])
          des2 <- rbind(des2, des1[p, ])
        }
      }
      dm2 <- dm2[2:nrow(dm2), 1]
      des2 <- des2[2:nrow(des2), ]

      dm2 <- diag(dm2)
      if (det(t(des2) %*% dm2 %*% des2) != 0) {
        a <- des2 %*% (solve(t(des2) %*% dm2 %*% des2)) %*% t(des2) %*% dm2
        acounter <- 1
        for (i in 1:ncodes) {
          for (j in 1:ncodes) {
            if (onezero[i, j] != 0) {
              zadjres[i, j] <- (freqs[i, j] - expfreq[i, j]) / sqrt(expfreq[i, j] *
                (1 - a[acounter, acounter]))
              acounter <- acounter + 1
            }
          }
        }
      } else {
        cat("\n\nWarning:\nA nonsingular matrix has been identified, which means that proper adjusted")
        cat("\nresiduals cannot be computed for this data, probably because there are no values for one")
        cat("\nor more codes. Try recoding using sequential integers, and redo the analyses.")
        cat("\nThe adjusted residuals that are printed below are based on equation 5 from")
        cat("\nBakeman & Quera (1995, p. 274), and are close approximations to the proper values.")
        cat("\nThe procedures recommended by Bakemen & Quera (1995, p. 276), Haberman (1979), and")
        cat("\nChristensen (1997) cannot be conducted with nonsingular matrices.\n\n")
      }

      for (i in 1:ncodes) {
        for (j in 1:ncodes) {
          if (onezero[i, j] == 0) {
            zadjres[i, j] <- 0
            yulesq[i, j] <- 0
            kappa[i, j] <- 0
            zkappa[i, j] <- 0
            pzadjres[i, j] <- 1
            pzkappa[i, j] <- 1
          }
        }
      }
    }
    # if (outfile == 1) {
      # outf <- freqs
    # } else if (outfile == 2) {
      # outf <- tprob
    # } else if (outfile == 3) {
      # outf <- zadjres
    # } else if (outfile == 4) {
      # outf <- yulesq
    # } else if (outfile == 5) {
      # outf <- kappa
    # }
#    out <- rbind(out, cbind(dindex, matrix(t(outf), 1, (ncodes * ncodes))))
    b <- as.matrix(labels[1:ncodes])
    if ((dindex < nrow(each)) && (output == "all")) {
      cat("\n\n\nGroup/Segment Number:", dindex, "\n")
    } else if (dindex == nrow(each)) {
      cat("\n\n\nPooled Data:\n\n")
      cat("\nNumber of cases/groups/segments:", nrow(each) - 1, "\n")
      
      dfsh <- ncodes * (ncodes - 1) * ((nrow(each) - 1) - 1)
      plrx2sh <- 1 - pchisq(abs(lrx2sh), dfsh)
      lchisq <- cbind(lrx2sh, dfsh, plrx2sh)
	  cat("\n\nLikelihood Ratio (Chi-Square) Test of Homogeneity/Stationarity = ",
	      round(lrx2sh,2),",  df = ",dfsh,",  p = ",round(plrx2sh,5),"\n\n",sep='')
    }

    if ((dindex == nrow(each)) || (output == "all")) {
      b <- as.matrix(labels[1:ncodes])
      bb <- rbind(b, "Totals")

      cfreq <- rbind(cbind(freqs, rowtots), cbind(coltots, ntrans))
      rownames(cfreq) <- bb
      colnames(cfreq) <- bb
      cat("\nCell Frequencies, Row & Column Totals, & N\n\n")
      print(cfreq)

      if ((adjacent == 0) || (adjacent == 2)) {
        cat("\n\nThe processed ONEZERO matrix appears below. In the ONEZERO matrix,")
        cat("\na 0 indicates a structural zero, and a 1 indicates that an expected ")
        cat("\ncell frequency will be estimated.\n")
        cat("\nONEZERO matrix:\n\n")
        rownames(onezero) <- b
        colnames(onezero) <- b
        print(onezero)
      }

      cat("\n\nExpected Values/Frequencies\n\n")
      rownames(expfreq) <- b
      colnames(expfreq) <- b
      print(round(expfreq, 2))

      cat("\n\nTransitional Probabilities\n\n")
      rownames(tprob) <- b
      colnames(tprob) <- b
      print(round(tprob, 2))

      if (adjacent == 1) {
        dft <- (ncodes - 1)^2
      } else {
        dft <- (ncodes - 1)^2 - (ncodes^2 - sum(onezero))
      }

	  plrx2t <- 1 - pchisq(abs(lrx2t), dft)
	  tlr <- cbind(lrx2t, dft, plrx2t)
	  cat("\n\nTablewise Likelihood Ratio (Chi-Square) test = ",round(lrx2t,2),
	      ",  df = ",dft,",  p = ",round(plrx2t,5),"\n",sep='')

      cat("\n\nAdjusted Residuals\n\n")
      rownames(zadjres) <- b
      colnames(zadjres) <- b
      print(round(zadjres, 2))

      cat("\n\nSignificance Levels for the Adjusted Residuals\n\n")
      rownames(pzadjres) <- b
      colnames(pzadjres) <- b
      print(round(pzadjres, 4))

      cat("\n\nYule's Q Values\n\n")
      rownames(yulesq) <- b
      colnames(yulesq) <- b
      print(round(yulesq, 2))

      cat("\n\nUnidirectional Kappas\n\n")
      rownames(kappa) <- b
      colnames(kappa) <- b
      print(round(kappa, 2))

      cat("\n\nz values for the Unidirectional Kappas\n\n")
      rownames(zkappa) <- b
      colnames(zkappa) <- b
      print(round(zkappa, 2))

      cat("\n\nSignificance Levels for the Unidirectional Kappas\n\n")
      rownames(pzkappa) <- b
      colnames(pzkappa) <- b
      print(round(pzkappa, 4))
    }
  }

  seqg_output <- list(
    freqs = freqs, expfreqs = expfreq, probs = tprob, chi = lchisq, adjres = zadjres,
    p = pzadjres, YulesQ = yulesq, kappas = kappa, z = zkappa, pk = pzkappa, output=output
  )

  return(invisible(seqg_output))
}
