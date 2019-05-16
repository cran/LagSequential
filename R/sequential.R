sequential <- function(data, labels = NULL, lag = 1, adjacent = TRUE,
                       onezero = NULL, tailed = 2, permtest = FALSE, nperms = 10) {

cat('\n\nLag Sequential Analysis\n')

  if (!adjacent) {
    adjacent <- 0
  } else if (adjacent && is.null(onezero)) {
    adjacent <- 1
  } else if (adjacent && !is.null(onezero)) {
    adjacent <- 2
  }

  if (!is.matrix(data)) data <- matrix(data, ncol = 1)


  # if data is a frequency transition matrix
  if (nrow(data) == ncol(data)) {
    datais <- 2
    ncodes <- ncol(data)
    freqs <- data
    if (is.null(labels)) {
      labels <- 1:ncodes
      for (lupe in 1:ncodes) labels[lupe] <- paste("Code", lupe)
    }
  }


  # if dataset is NOT a frequency transition matrix
  if (!(nrow(data) == ncol(data))) {
    datais <- 1

    data <- matrix(data, ncol = 1)

    # are all data values numeric? any problems?
    if ((all(sapply(data, is.numeric))) == TRUE) {
      codesmin <- min(data)
      codesmax <- max(data)
      codefreqs <- table(data)
      cat("\n\nThe code frequencies:\n\n")
      print(codefreqs)
      if (codesmin != 1 || (codesmax > length(codefreqs)) ||
          (min(codefreqs) == 0)) {
        warning("The entered data is numeric but there is a problem:",
                "\n\t-- the minimum code value should be 1,",
                "\n\t-- the set of possible code values should be consecutive",
                "integers",
                "\n\t-- all code frequencies should be 1",
                "\nAt least one of these conditions has not been met",
                " which will cause problems.")
      }

      ncodes <- max(data)

      if (is.null(labels)) {
        labels <- 1:ncodes
        for (lupe in 1:ncodes) labels[lupe] <- paste(" Code", lupe)
      }
    }

    # if any data values are characters, treat them all as strings & provide numeric values for the analyses
    if (any(sapply(data, is.character))) {
      labels <- unique(data)
      for (lupe in 1:length(data)) data[lupe, 1] <- which(labels == data[
          lupe,
          1
        ], arr.ind = F)
      data <- as.matrix(as.numeric(data))
      ncodes <- max(data)
    }

    # transitional frequency matrix.
    freqs <- matrix(0, ncodes, ncodes)
    for (c in 1:nrow(data)) {
      if (c + lag <= nrow(data)) {
        freqs[data[c], data[c + lag]] <- freqs[data[c], data[c + lag]] + 1
      }
    }
  }

  # Warning message for specification error.
  if ((adjacent == 0) & (any(diag(freqs) == 0))) {
    warning("\n\nYou have indicated that adjacent codes can never repeat", 
            "\n(ajdacent = 0), yet repeating values have been found in the data.",
            "\nSee the main diagonal frequency matrix.",
            "\nThis will result in faulty computations for LRX2,",
            "\nz values, and adjusted residuals.\n\n")
  }

  # initializing.
  lrx2t <- matrix(0)
  rowtots <- matrix(rowSums(freqs), ncol = 1)
  coltots <- matrix(colSums(freqs), nrow = 1)
  ntrans <- sum(rowtots)
  prows <- rowtots / ntrans
  pcols <- coltots / ntrans
  tprob <- matrix(-9999, ncodes, ncodes)
  et <- matrix(-9999, ncodes, ncodes)
  expfreq <- matrix(-9999, ncodes, ncodes)
  zadjres <- matrix(-9999, ncodes, ncodes)
  pzadjres <- matrix(1, ncodes, ncodes)
  yulesq <- matrix(-9999, ncodes, ncodes)
  var <- matrix(-9999, ncodes, ncodes)
  min <- matrix(-9999, ncodes, ncodes)
  kappa <- matrix(-9999, ncodes, ncodes)
  zkappa <- matrix(-9999, ncodes, ncodes)
  pzkappa <- matrix(1, ncodes, ncodes)
  signs <- matrix(0, ncodes, ncodes)
  n <- ntrans + 1
  nr <- rowtots

  if (datais == 1) {
    nr[data[nrow(data), 1]] <- nr[data[nrow(data), 1]] + 1
  }


  for (i in 1:ncodes) {
    for (j in 1:ncodes) {

      # Note: more refined computations for when adjacent codes cannot repeat appear below,
      # after the above 2 loops are completed.
      if (adjacent == 0 & (ntrans - rowtots[i]) > 0) {
        pcols[j] <- coltots[j] / (ntrans - rowtots[i])
      }
      if (adjacent == 0 & (ntrans - rowtots[j]) > 0) {
        expfreq[i, j] <- (rowtots[i] * coltots[j]) / (ntrans - rowtots[j])
      }
      if (adjacent == 0 & (n - nr[i]) > 0) {
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
      if (freqs[i, j] > 0 & expfreq[i, j] > 0) {
        lrx2t <- lrx2t + 2 * (freqs[i, j] * log(freqs[i, j] / expfreq[i, j]))
      }

      # adjusted residuals (z values) & sig levels
      if ((expfreq[i, j] * (1 - pcols[j]) * (1 - prows[i])) > 0) {
        zadjres[i, j] <- (freqs[i, j] - expfreq[i, j]) / sqrt(expfreq[i, j] *
          (1 - pcols[j]) * (1 - prows[i]))
        pzadjres[i, j] <- (1 - pnorm(abs(zadjres[i, j]))) * tailed
      }

      # Yule's Q.
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

      # signs
      if (freqs[i, j] > expfreq[i, j]) {
        signs[i, j] <- 1
      } else if (freqs[i, j] < expfreq[i, j]) {
        signs[i, j] <- (-1)
      }
    }
  }

  if ((adjacent == 0) || (adjacent == 2)) {

    # maximum likelihood estimation of the expected cell frequencies using iterative proportional fitting (Wickens, 1989, pp. 107-112).

    rsumsf <- rowSums(freqs)
    csumsf <- colSums(freqs)

    if (is.null(onezero)) {
      onezero <- matrix(1, ncodes, ncodes)
      diag(onezero) <- 0
    }

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

    cat("\nMaximum likelihood estimation of the expected cell frequencies using iterative proportional fitting ")
    if ((max(rdiffs) < 1e-04) & (max(cdiffs) < 1e-04)) {
      cat(
        "converged after the following number of iterations:", ipfloop,
        "\n\n"
      )
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
          lrx2t <- lrx2t + 2 * (freqs[i, j] * log(freqs[i, j] / expfreq[i, j]))
        }
      }
    }

    # adjusted residuals for matrices with structural zeros (Christensen, 1997, p. 357).

    # constructing the design matrix.
    x <- matrix(1, ncodes^2, 1)
    y <- matrix(0, ncodes^2, ncodes - 1)
    z <- matrix(0, ncodes^2, ncodes - 1)
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
    for (pp in 1:(ncodes^2)) {
      if (onezero2[pp] == 1) {
        dm2 <- rbind(dm2, dm1[pp])
        des2 <- rbind(des2, des1[pp, ])
      }
    }
    dm2 <- dm2[2:nrow(dm2), 1]
    des2 <- des2[2:nrow(des2), ]

    dm2 <- diag(dm2)
    if (det(t(des2) %*% dm2 %*% des2) != 0) {
      zadjres <- matrix(0, ncodes, ncodes)
      a <- des2 %*% (solve(t(des2) %*% dm2 %*% des2)) %*% t(des2) %*% dm2
      acounter <- 1
      for (i in 1:ncodes) {
        for (j in 1:ncodes) {
          if (onezero[i, j] != 0) {
            zadjres[i, j] <- (freqs[i, j] - expfreq[i, j]) / sqrt(expfreq[
              i,
              j
            ] * (1 - a[acounter, acounter]))
            acounter <- acounter + 1
          }
        }
      }
    } else {
      warning("\n\nA nonsingular matrix has been identified, which means that proper",
          "\nadjusted residuals cannot be computed for this data, probably",
          "\nbecause there are no values for one or more codes. Try recoding",
          "\nusing sequential integers, and redo the analyses. The adjusted",
          "\nresiduals that are printed below are based on equation 5 from",
          "\nBakemand & Quera (1995, p. 274), and are close approximations",
          "\nto the proper values. The procedures recommended by Bakemen &",
          "\nQuera (1995, p. 276), Haberman (1979), and Christensen (1997)",
          "\ncannot be conducted with nonsingular matrices.\n\n")
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

  b <- labels[1:ncodes]
  bb <- c(b, "Totals")
  cat("\n\nRequested 'tail' (1 or 2) for Significance Tests =", tailed, "\n")

  cfreq <- rbind(cbind(freqs, rowtots), cbind(coltots, ntrans))
  rownames(cfreq) <- bb
  colnames(cfreq) <- bb
  cat("\n\nCell Frequencies, Row & Column Totals, & N\n\n")
  print(cfreq)

  if (adjacent == 0 || adjacent == 2) {
    cat("\nThe processed ONEZERO matrix appears below. In the ONEZERO matrix,",
        "\na 0 indicates a structural zero, and a 1 indicates that an expected cell",
        "\nfrequency will be estimated.",
        "\n\nONEZERO matrix:\n")
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
    df <- (ncodes - 1)^2
  } else {
    df <- (ncodes - 1)^2 - (ncodes^2 - sum(onezero))
  }

  plrx2t <- 1 - pchisq(abs(lrx2t), df)
  tlr <- cbind(lrx2t, df, plrx2t)
  cat("\n\nTablewise Likelihood Ratio (Chi-Square) test = ",round(lrx2t,2),
      ",  df = ",df,",  p = ",round(plrx2t,5),"\n",sep='')

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




  # Permutation tests of significance

  if (permtest && datais == 2) {
    warning("\n\nYou have requested permutation tests of significance be computed",
            "\n(permtest = TRUE) but the data is a frequency transition matrix.",
            "\nPermutation tests can only be computed when data is a sequence",
            "\nof codes. Permutation tests will not be performed.\n\n")
  }

  if (permtest && datais == 1) {
    obs2 <- matrix(t(freqs), 1, nrow(freqs) * ncol(freqs))
    signs2 <- matrix(t(signs), 1, nrow(freqs) * ncol(freqs))
    sigs <- matrix(1, 1, nrow(freqs) * ncol(freqs))

      results <- matrix(-9999, nperms, nrow(freqs) * ncol(freqs))

      for (perm in 1:nperms) {

        # permuting the sequences; algorithm from Castellan 1992.

        # when adjacent codes may be the same.
        datap <- data
        if (adjacent == 1) {
          for (i in 1:(nrow(datap) - 1)) {
            kay <- as.integer((nrow(datap) - i + 1) * runif(1) + 1) + i - 1
            d <- datap[i]
            datap[i] <- datap[kay]
            datap[kay] <- d
          }
          datap <- matrix(sample(data), nrow(data), 1)
        }

        # when adjacent codes may NOT be the same.
        if (adjacent == 0) {
          datap <- rbind(0, data, 0)
          for (i in 2:(nrow(datap) - 2)) {
            limit <- 10000
            for (j in 1:limit) {
              kay <- as.integer(((length(datap) - 1) - i + 1) * runif(1) + 1) + i - 1
              if ((datap[i - 1] != datap[kay]) & (datap[i + 1] != datap[kay]) &
                 (datap[kay - 1] != datap[i]) & (datap[kay + 1] != datap[i])) {
                break
              }
            }
            d <- datap[i]
            datap[i] <- datap[kay]
            datap[kay] <- d
          }
          datap <- matrix(datap[2:(nrow(datap) - 1), ], ncol = 1)
        }

        # transitional frequency matrix for permuted data.
        freqsp <- matrix(0, ncodes, ncodes)
        for (c in 1:nrow(datap)) {
          if (c + lag <= nrow(datap)) {
            freqsp[datap[c], datap[c + lag]] <- freqsp[datap[c], datap[c + lag]] + 1
          }
        }
        results[perm, ] <- matrix(t(freqsp), 1, nrow(freqs) * ncol(freqs))
      }


      # one-tailed
      for (j in 1:ncol(results)) {
	  counter <- 0
          for (i in 1:nrow(results)) {
            if ((results[i, j] >= obs2[j]) & (signs2[j] > 0)) { 
            	counter <- counter + 1
            } else if ((results[i, j] <= obs2[j]) & (signs2[j] < 0)) {
            	counter <- counter + 1
            }
          }
          if (signs2[j] != 0) sigs[1, j] <- counter / nperms
        }
        

    cat("\n\nData Permutation Significance Levels (number of permutations = ", nperms, ")\n\n",sep='')
    sigs <- t(matrix(sigs, ncodes, ncodes))
    rownames(sigs) <- b
    colnames(sigs) <- b    
    print(sigs)

  }

  seq_output <- list(
    freqs = freqs, expfreqs = expfreq, probs = tprob, chi = tlr,
    adjres = zadjres, p = pzadjres, YulesQ = yulesq, kappas = kappa, z = zkappa,
    pk = pzkappa
  )

  return(invisible(seq_output))
}
