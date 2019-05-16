# BIDIRECTIONAL

# This function tests the bidirectional dependence of behaviors i to j,
# and j to i, an additive sequential pattern described by Wampold and
# Margolin (1982) and Wampold (1989, 1992).


bidirectional <- function(data, labels = NULL, lag = 1, adjacent = TRUE, tailed = 1,
                          permtest = FALSE, nperms = 10) {

  cat('\n\nLag Sequential Analysis Tests for Bidirectional Dependence\n')
  
    if (is.matrix(data) == FALSE) data <- matrix(data, ncol = 1)


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


  # if data is NOT a frequency transition matrix
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
      if (codesmin != 1 | (codesmax > length(codefreqs)) | (min(codefreqs) == 0)) {
        cat("\n\nThe entered data is numeric, but there is a problem:")
        cat("\n    -- the minimum code value should be 1,")
        cat("\n    -- the set of possible code values should be consecutive integers, &")
        cat("\n    -- all code frequencies should > 1")
        cat("\nAt least one of these conditions has not been met, which will cause problems.")
      }

      ncodes <- max(data)

      if (is.null(labels)) {
        labels <- 1:ncodes
        for (lupe in 1:ncodes) labels[lupe] <- paste("Code", lupe)
      }
    }

    # if any data values are characters, treat them all as strings & provide numeric values for the analyses
    if ((any(sapply(data, is.character))) == TRUE) {
      labels <- unique(data)
      for (lupe in 1:length(data)) data[lupe, 1] <- which(labels == data[lupe, 1], arr.ind = F)
      data <- as.matrix(as.numeric(data))
      ncodes <- max(data)
    }

    # transitional frequency matrix.
    freqs <- matrix(0, ncodes, ncodes)
    for (c in 1:nrow(data)) {
      if (c + lag <= nrow(data)) freqs[data[c], data[c + lag]] <- freqs[data[c], data[c + lag]] + 1
    }
  }

  # initializing
  ett <- matrix(-9999, ncodes, ncodes)
  var <- matrix(-9999, ncodes, ncodes)
  min <- matrix(-9999, ncodes, ncodes)
  kappa <- matrix(-9999, ncodes, ncodes)
  zkappa <- matrix(-9999, ncodes, ncodes)
  pzkappa <- matrix(1, ncodes, ncodes)
  signs <- matrix(0, ncodes, ncodes)
  obs <- matrix(0, ncodes, ncodes)
  rowtots <- matrix(rowSums(freqs))
  coltots <- matrix(colSums(freqs), ncol = ncodes)
  ntrans <- sum(rowtots)
  n <- ntrans + 1
  nr <- rowtots

  if (datais == 1) nr[data[nrow(data), 1]] <- nr[data[nrow(data), 1]] + 1


  for (i in 1:ncodes) {
    for (j in 1:ncodes) {
      if ((nr[i] > 0) & (nr[j] > 0)) {
        obs[i, j] <- freqs[i, j] + freqs[j, i]
        ett[i, j] <- 2 * nr[i] * nr[j] / n
        var[i, j] <- 2 * nr[i] * nr[j] * (nr[i] * nr[j] + (n - nr[i]) * (n - nr[j]) -
          n) / (n^2 * (n - 1))
        zkappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - ett[i, j]) / sqrt(var[
          i,
          j
        ])
        pzkappa[i, j] <- (1 - pnorm(abs(zkappa[i, j]))) * tailed
        if (nr[i] <= nr[j]) {
          min[i, j] <- nr[i]
        } else {
          min[i, j] <- nr[j]
        }
        kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - ett[i, j]) / (2 * min[i, j] -
          ett[i, j])
        if (kappa[i, j] < 0) {
          kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - ett[i, j]) / ett[i, j]
        }
        if (nr[i] == nr[j]) {
          kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - ett[i, j]) / ((2 * nr[j] -
            1) - ett[i, j])
        }
        # signs.
        if (kappa[i, j] > 0) {
          signs[i, j] <- 1
        } else if (kappa[i, j] < 0) {
          signs[i, j] <- -1
        }
      }
    }
  }

  b <- labels[1:ncodes]
  bb <- c(b, "Totals")
  cfreqs <- rbind(cbind(freqs, rowtots), cbind(coltots, sum(rowtots)))
  rownames(cfreqs) <- bb
  colnames(cfreqs) <- bb
  cat("\n\nCell Frequencies, Row & Column Totals, & N\n\n")
  print(cfreqs)

  rownames(obs) <- b
  colnames(obs) <- b
  cat("\n\nObserved Bidirectional Frequencies\n\n")
  print(obs)

  rownames(ett) <- b
  colnames(ett) <- b
  cat("\n\nExpected Bidirectional Frequencies\n\n")
  print(round(ett, 2))

  rownames(kappa) <- b
  colnames(kappa) <- b
  cat("\n\nBidirectional Kappas\n\n")
  print(round(kappa, 2))

  rownames(zkappa) <- b
  colnames(zkappa) <- b
  cat("\n\nz values for the bidirectional Kappas\n\n")
  print(round(zkappa, 2))

  cat("\n\nRequested 'tail' (1 or 2) for Significance Tests =", tailed, "\n")

  rownames(pzkappa) <- b
  colnames(pzkappa) <- b
  cat("\n\nSignificance Levels for the Bidirectional Kappas\n\n")
  print(round(pzkappa, 4))



  # Permutation tests of significance
  if (permtest && datais == 1) {
    obs2 <- matrix(t(obs), 1, (nrow(freqs) * ncol(freqs)))
    obs22 <- matrix(t((ett - (obs - ett))), 1, (nrow(freqs) * ncol(freqs)))
    signs2 <- matrix(t(signs), 1, (nrow(freqs) * ncol(freqs)))
    sigs <- matrix(1, 1, (nrow(freqs) * ncol(freqs)))

      results <- matrix(-9999, nperms, (nrow(freqs) * ncol(freqs)))

      for (perm in 1:nperms) {

        # permuting the sequences; algorithm from Castellan 1992.

        # when adjacent codes may be the same.
        datap <- data
        if (adjacent) {
          for (i in 1:(nrow(datap) - 1)) {
            kay <- as.integer((nrow(datap) - i + 1) * runif(1) + 1) + i - 1
            d <- datap[i]
            datap[i] <- datap[kay]
            datap[kay] <- d
          }
        }

        # when adjacent codes may NOT be the same.
        if (!adjacent) {
          datap <- rbind(0, data, 0)
          for (i in 2:(nrow(datap) - 2)) {
            limit <- 10000
            for (j in 1:limit) {
              kay <- as.integer(((nrow(datap) - 1) - i + 1) * runif(1) + 1) + i - 1
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

        # bidirectional frequency matrix for permuted data.
        obsp <- matrix(0, ncodes, ncodes)
        for (i in 1:ncodes) {
          for (j in 1:ncodes) {
            obsp[i, j] <- freqsp[i, j] + freqsp[j, i]
          }
        }

        results[perm, ] <- matrix(t(obsp), 1, nrow(freqs) * ncol(freqs))
      }

      # one-tailed.
      if (tailed == 1) {
        for (j in 1:ncol(results)) {
          counter <- 0
          for (i in 1:nrow(results)) {
            if ((results[i, j] >= obs2[j]) & (signs2[j] > 0)) {
              counter <- counter + 1
            } else if ((results[i, j] <= obs2[j]) & (signs2[j] < 0)) {
              counter <- counter + 1
            }
          }
          if (signs2[j] != 0) {
            sigs[1, j] <- counter / nperms
          }
        }
      }

    cat("\n\nData Permutation Significance Levels (number of permutations = ", nperms,")\n\n",sep='')
    sigs <- t(matrix(sigs, ncodes, ncodes))
    rownames(sigs) <- b
    colnames(sigs) <- b    
    print(sigs)
  }

  bid_output <- list(freqs = freqs, bifreqs = obs, expbifreqs = ett, kappas = kappa, z = zkappa, pk = pzkappa)

  return(invisible(bid_output))
}
