# TWOCELLS

# This function simultaneously tests the unidirectional dependence of i
# to j, and the unidirectional dependence of k to L, an additive
# pattern described by Wampold and Margolin (1982) and Wampold (1989,
# 1992).

twocells <- function(data, i, j, k, L, labels = NULL, lag = 1,
                     adjacent = TRUE, tailed = 1, permtest = FALSE, nperms = 10) {
                     	
cat('\n\nLag Sequential Analysis "Two Cells" Tests\n')
cat('\n\ne.g., of the unidirectional dependence of i to j and the unidirectional dependence of k to L\n')


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
  if ((nrow(data) == ncol(data)) == FALSE) {
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


  # when the supplied values for i, j, k, or L are not integers i.e., are strings
  if ((any(sapply(c(i, j, k, L), is.character))) == TRUE) {
    i <- which(labels == i, arr.ind = F)
    j <- which(labels == j, arr.ind = F)
    k <- which(labels == k, arr.ind = F)
    L <- which(labels == L, arr.ind = F)
  }


  rowtots <- matrix(rowSums(freqs))
  coltots <- matrix(colSums(freqs), ncol = ncodes)
  ntrans <- sum(rowtots)
  n <- ntrans + 1
  nr <- rowtots

  if (datais == 1) nr[data[nrow(data), 1]] <- nr[data[nrow(data), 1]] + 1


  nr[data[n]] <- nr[data[n]] + 1
  ett <- (nr[i] * nr[j] + nr[k] * nr[L]) / n
  var <- (nr[i] * nr[j] * (n - nr[i]) * (n - nr[j]) + nr[k] * nr[L] * (n - nr[k]) *
    (n - nr[L]) + 2 * nr[i] * nr[j] * nr[k] * nr[L]) / (n^2 * (n - 1))
  zkappa <- ((freqs[i, j] + freqs[k, L]) - ett) / sqrt(var)
  pzkappa <- (1 - pnorm(abs(zkappa))) * tailed
  if (nr[i] <= nr[j]) {
    minij <- nr[i]
  } else {
    minij <- nr[j]
  }
  if (nr[k] <= nr[L]) {
    minkL <- nr[k]
  } else {
    minkL <- nr[L]
  }
  kappa <- ((freqs[i, j] + freqs[k, L]) - (nr[i] * nr[j] + nr[k] * nr[L]) / n) / (minij +
    minkL - ((nr[i] * nr[j] + nr[k] * nr[L]) / n))
  if (kappa < 0) {
    kappa <- ((freqs[i, j] + freqs[k, L]) - (nr[i] * nr[j] + nr[k] * nr[L]) / n) / (((nr[i] *
      nr[j] + nr[k] * nr[L]) / n))
  }
  b <- labels[1:ncodes]
  bb <- c(b, "Totals")

  cfreqs <- rbind(cbind(freqs, rowtots), cbind(coltots, sum(rowtots)))
  rownames(cfreqs) <- bb
  colnames(cfreqs) <- bb
  cat("\n\nCell Frequencies, Row & Column Totals, & N\n\n")
  print(cfreqs)

  cat("\n\nSimultaneous Two-Cell Test for the following code values:
	   \n     Code i =",labels[i],
	  "\n\n     Code j =",labels[j],
	  "\n\n     Code k =",labels[k],
	  "\n\n     Code L =",labels[L])

#  cat("\n\nRequested 'tail' (1 or 2) for Significance Tests =", tailed, "\n")

  cat("\n\n\nkappa =",round(kappa,2),"   z =",round((zkappa),3),"   p =",round(pzkappa,5),"\n\n")




  # Permutation tests of significance

  if (permtest && datais == 1) {
  	
    obs2 <- freqs[i, j] + freqs[k, L]
    obs22 <- ett - (obs2 - ett)
    if (kappa > 0) {
      sign <- 1
    } else if (kappa < 0) {
      sign <- -1
    } else {
      sign <- 0
    }
    sigs <- matrix(1, 1, 1)

      results <- matrix(-9999, nperms, 1)

      for (perm in 1:nperms) {

        # permuting the sequences; algorithm from Castellan 1992.

        # when adjacent codes may be the same.
        datap <- data
        if (adjacent) {
          for (ii in 1:(nrow(datap) - 1)) {
            kay <- as.integer((nrow(datap) - ii + 1) * runif(1) + 1) + ii - 1
            d <- datap[ii]
            datap[ii] <- datap[kay]
            datap[kay] <- d
          }
        }

        # when adjacent codes may NOT be the same.
        if (!adjacent) {
          datap <- rbind(0, data, 0)
          for (ii in 2:(nrow(datap) - 2)) {
            limit <- 10000
            for (jj in 1:limit) {
              kay <- as.integer(((nrow(datap) - 1) - ii + 1) * runif(1) + 1) + ii - 1
              if ((datap[ii - 1] != datap[kay]) & (datap[ii + 1] != datap[kay]) &
                (datap[kay - 1] != datap[ii]) & (datap[kay + 1] != datap[ii])) {
                break
              }
            }
            d <- datap[ii]
            datap[ii] <- datap[kay]
            datap[kay] <- d
          }
          datap <- matrix(datap[2:(nrow(datap) - 1), ], ncol = 1)
        }

        # transitional frequency matrix for permuted data
        freqsp <- matrix(0, ncodes, ncodes)
        for (c in 1:nrow(datap)) {
          if (c + lag <= nrow(datap)) {
            freqsp[datap[c], datap[c + lag]] <- freqsp[datap[c], datap[c + lag]] + 1
          }
        }

        # two-cell frequency for permuted data.
        obsp <- freqsp[i, j] + freqsp[k, L]

        results[perm, 1] <- obsp
      }

      # one-tailed.
      if (tailed == 1) {
        counter <- 0
        for (ii in 1:nrow(results)) {
          if (results[ii] >= obs2 & sign > 0) {
            counter <- counter + 1
          } else if (results[ii] <= obs2 & sign < 0) {
            counter <- counter + 1
          }
        }
        if (sign != 0) {
          sigs[1, 1] <- counter / nperms
        }
      }

    cat("\nData Permutation Significance Level (for ", nperms, " permutations) = ", sigs, "\n\n\n", sep='')

  }

  twocell_output <- list(
    freqs = freqs, twocellfreq = (freqs[i, j] + freqs[k, L]), expfreqs = ett,
    kappa = kappa, z = zkappa, pk = pzkappa
  )

  return(invisible(twocell_output))
}