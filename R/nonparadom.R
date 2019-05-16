# NONPARADOM

# This function tests for nonparallel dominance, a form of asymmetry in
# predictability, between i to j, and k to L, as described by Wampold
# (1984, 1989, 1992).

nonparadom <- function(data, i, j, k, L, labels = NULL, lag = 1, adjacent = TRUE, 
                       tailed = 1, permtest = FALSE, nperms = 10) {

cat('\n\nLag Sequential Analysis Tests for Nonparallel Dominance\n')
 
 if (is.matrix(data) == FALSE) data <- matrix(data, ncol = 1)


# if data is a frequency transition matrix
if (nrow(data)==ncol(data)) {
	datais <- 2
	ncodes <- ncol(data)
	freqs <- data	
	if (is.null(labels)) {
		 labels <- 1:ncodes		 
		 for (lupe in 1:ncodes)  labels[lupe] <- paste('Code',lupe)
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
		cat('\n\nThe code frequencies:\n\n'); print(codefreqs)		
			if ( codesmin !=1 | (codesmax > length(codefreqs)) | (min(codefreqs) == 0)) {
				cat('\n\nThe entered data is numeric, but there is a problem:')
				cat('\n    -- the minimum code value should be 1,')
				cat('\n    -- the set of possible code values should be consecutive integers, &')
				cat('\n    -- all code frequencies should > 1')
				cat('\nAt least one of these conditions has not been met, which will cause problems.')
			}
			
		ncodes <- max(data)

		if (is.null(labels)) {
			labels <- 1:ncodes		 
			for (lupe in 1:ncodes)  labels[lupe] <- paste('Code',lupe)
		}	
	}

	# if any data values are characters, treat them all as strings & provide numeric values for the analyses
	if ((any(sapply(data, is.character))) == TRUE) {		
		labels <- unique(data)
		for (lupe in 1:length(data))  data[lupe,1] <- which(labels == data[lupe,1], arr.ind = F)
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
if ((any(sapply(c(i,j,k,L), is.character))) == TRUE) {		
	i <- which(labels == i, arr.ind = F)
	j <- which(labels == j, arr.ind = F)
	k <- which(labels == k, arr.ind = F)
	L <- which(labels == L, arr.ind = F)
}


	pd <- matrix(-9999, ncodes, ncodes)
	et <- matrix(-9999, ncodes, ncodes)
	rowtots <- matrix(rowSums(freqs))
	coltots <- matrix(colSums(freqs), ncol = ncodes)
    ntrans <- sum(rowtots)
    n <- ntrans + 1
    nr <- rowtots
  
    if (datais == 1) nr[data[nrow(data), 1]] <- nr[data[nrow(data), 1]] + 1

	prow <- nr/sum(nr)

	for (iindex in 1:ncodes) {
		for (jindex in 1:ncodes) {
			if (nr[iindex] > 0 & nr[jindex] > 0 & prow[jindex] > 0) {
				pd[iindex, jindex] <- ((freqs[iindex, jindex]/nr[iindex]) - prow[jindex])/prow[jindex]
				if (nr[iindex] > 0 & nr[jindex] > 0) {
					et[iindex, jindex] <- (nr[iindex] * nr[jindex])/n
				}
			}
		}
	}

	# 96

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

	# kappa.
	kappa <- -9999
	case <- -9999
	if (freqs[i, j] == et[i, j]) {
		kappa <- 0
		case <- 0
	}

	if (nr[i] > 0 & nr[j] > 0 & nr[k] > 0 & nr[L] > 0) {

		# Wampold's 1st case.
		if (freqs[i, j] > et[i, j] & freqs[k, L] >= et[k, L]) {
			kappa <- (nr[k] * nr[L] * freqs[i, j] - nr[i] * nr[j] * freqs[k, L])/(nr[k] *
				nr[L] * minij - (nr[i] * nr[j] * nr[k] * nr[L]/n))
			if (kappa < 0) {
				kappa <- (nr[k] * nr[L] * freqs[i, j] - nr[i] * nr[j] * freqs[k, L])/(nr[i] *
					nr[j] * minkL - (nr[i] * nr[j] * nr[k] * nr[L]/n))
			}
			case <- 1
		}

		# Wampold's 2nd case.
		if (freqs[i, j] < et[i, j] & freqs[k, L] <= et[k, L]) {
			kappa <- (nr[k] * nr[L] * freqs[i, j] - nr[i] * nr[j] * freqs[k, L])/((-1) *
				(nr[i] * nr[j] * nr[k] * nr[L]/n))
			case <- 2
		}

		# Wampold's 3rd case.
		if (freqs[i, j] > et[i, j] & freqs[k, L] <= et[k, L]) {
			kappa <- (nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L] - 2 *
				(nr[i] * nr[j] * nr[k] * nr[L]/n))/(nr[k] * nr[L] * minij - (nr[i] * nr[j] *
				nr[k] * nr[L]/n))
			if (kappa == 0) {
				kappa <- (nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L] -
					2 * (nr[i] * nr[j] * nr[k] * nr[L]/n))/(nr[i] * nr[j] * nr[k] * nr[L]/n)
			}
			case <- 3
		}

		# Wampold's 4th case.
		if (freqs[i, j] < et[i, j] & freqs[k, L] >= et[k, L]) {
			kappa <- (nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L] - 2 *
				(nr[i] * nr[j] * nr[k] * nr[L]/n))/((-1) * (nr[i] * nr[j] * nr[k] * nr[L]/n))
			if (kappa < 0) {
				kappa <- (nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L] -
					2 * (nr[i] * nr[j] * nr[k] * nr[L]/n))/((nr[i] * nr[j] * nr[k] * nr[L]/n) -
					nr[i] * nr[j] * minkL)
			}
			case <- 4
		}
	}

	# observed frequency, expected frequency, variance, & z.
	zeqk <- 9999
	obs <- -9999
	ett <- -9999
	zkappa <- -9999
	pzkappa <- -9999
	if (pd[i, j] != -9999 & pd[k, L] != -9999) {

		# same direction.
		if ((pd[i, j] >= 0 & pd[k, L] >= 0) || (pd[i, j] <= 0 & pd[k, L] <= 0)) {
			obs <- nr[k] * nr[L] * freqs[i, j] - nr[i] * nr[j] * freqs[k, L]
			ett <- 0
			var <- (nr[i] * nr[j] * nr[k] * nr[L] * (n * nr[i] * nr[j] + n * nr[k] *
				nr[L] - nr[i] * nr[j] * nr[k] - nr[i] * nr[j] * nr[L] - nr[i] * nr[k] *
				nr[L] - nr[j] * nr[k] * nr[L]))/(n * (n - 1))
			zkappa <- obs/sqrt(var)
			zeqk <- zkappa
			# different directions
			} else if ((pd[i, j] <= 0 & pd[k, L] >= 0) || (pd[i, j] >= 0 & pd[k, L] <= 0)) {
			obs <- nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L]
			ett <- 2 * nr[i] * nr[j] * nr[k] * nr[L]/n
			var <- (nr[i] * nr[j] * nr[k] * nr[L] * (4 * nr[i] * nr[j] * nr[k] * nr[L] +
				n^2 * nr[i] * nr[j] + n^2 * nr[k] * nr[L] - n * nr[i] * nr[j] * nr[k] -
				n * nr[i] * nr[j] * nr[L] - n * nr[i] * nr[k] * nr[L] - n * nr[j] * nr[k] *
				nr[L]))/(n^2 * (n - 1))
			zkappa <- (obs - ett)/sqrt(var)
			if (((pd[i, j] == 0 || pd[j, i] == 0)) & zeqk < zkappa) {
				zkappa <- zeqk
			}
		}
		pzkappa <- (1 - pnorm(abs(zkappa))) * tailed
	}

	if (kappa > 0) {
		sign <- 1
	} else if (kappa < 0 & kappa != -9999) {
		sign <- (-1)
	} else {
		sign <- 0
	}

	b <- labels[1:ncodes]
	bb <- c(b, "Totals")

	cfreqtotn <- rbind(cbind(freqs, rowtots), cbind(coltots, sum(rowtots)))
	rownames(cfreqtotn) <- bb
	colnames(cfreqtotn) <- bb
	cat("\n\nCell Frequencies, Row & Column Totals, & N\n\n")
	print(cfreqtotn)

	rownames(et) <- b
	colnames(et) <- b
	cat("\n\nExpected Values/Frequencies\n\n")
	print(round(et,2))

	cat("\n\nNonparallel Dominance test for the following code values:
		   \n     Code i =",labels[i],
		"\n\n     Code j =",labels[j],
		"\n\n     Code k =",labels[k],
		"\n\n     Code L =",labels[L])

	cat("\n\n\nSequential Dominance 'Case' type (Wampold, 1989):\n")
	if (case == 1) cat("\n     Case 1:",labels[i],"increases",labels[j],"and",labels[k],"increases",labels[L])
	if (case == 2) cat("\n     Case 2:",labels[i],"decreases",labels[j],"and",labels[k],"decreases",labels[L])
	if (case == 3) cat("\n     Case 3:",labels[i],"increases",labels[j],"and",labels[k],"decreases",labels[L])
	if (case == 4) cat("\n     Case 4:",labels[i],"decreases",labels[j],"and",labels[k],"increases",labels[L])


#	cat("\n\nRequested 'tail' (1 or 2) for Significance Tests =",tailed,"\n\n")

	cat("\n\n\nkappa =",round(kappa,2),"   z =",round((zkappa*sign),3),"   p =",round(pzkappa,5),"\n\n")


	# Permutation tests of significance.
	if (permtest && datais == 1 && obs != -9999 && ett != -9999 & case > 0) {

		obs2 <- obs
		obs22 <- ett - (obs2 - ett)
		sigs <- matrix(1, 1, 1)

			results <- matrix(-9999, nperms, 1)

			for (perm in 1:nperms) {

				# permuting the sequences; algorithm from Castellan 1992.

				# when adjacent codes may be the same.
				datap <- data
				if (adjacent) {
					for (iindex in 1:(nrow(datap) - 1)) {
						kay <- as.integer((nrow(datap) - iindex + 1) * runif(1) + 1) + iindex - 1
						d <- datap[iindex]
						datap[iindex] <- datap[kay]
						datap[kay] <- d
					}
				}

				# when adjacent codes may NOT be the same.
				if (!adjacent) {
					datap <- rbind(0, data, 0)
					for (iindex in 2:(nrow(datap) - 2)) {
						limit <- 10000
						for (jindex in 1:limit) {
							kay <- as.integer(((nrow(datap) - 1) - iindex + 1) * runif(1) + 1) + iindex - 1
							if ((datap[iindex - 1] != datap[kay]) & (datap[iindex + 1] != datap[kay]) &
								(datap[kay - 1] != datap[iindex]) & (datap[kay + 1] != datap[iindex])) {
								break
							}
						}
						d <- datap[iindex]
						datap[iindex] <- datap[kay]
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

				# nonparallel dominance frequency for permuted data
				np <- nrow(datap)
				nrp <- rowSums(freqsp)
				nrp[datap[np]] <- nrp[datap[np]] + 1
				if (case == 1 || case == 2) {
					obsp <- nrp[k] * nrp[L] * freqsp[i, j] - nrp[i] * nrp[j] * freqsp[k,L]
				} else if (case == 3 || case == 4) {
					obsp <- nrp[k] * nrp[L] * freqsp[i, j] + nrp[i] * nrp[j] * freqsp[k,L]
				}

				results[perm, 1] <- obsp
			}


			# one-tailed.
				counter <- 0
				for (iindex in 1:nrow(results)) {
					if (case == 1 || case == 3) {
						if (sign > 0 & results[iindex] >= obs2) {
							counter <- counter + 1
						} else if (sign < 0 & results[iindex] <= obs2) {
							counter <- counter + 1
						}
					}
					if (case == 2 || case == 4) {
						if (sign > 0 & results[iindex] <= obs2) {
							counter <- counter + 1
						} else if (sign < 0 & results[i] >= obs2) {
							counter <- counter + 1
						}
					}
				}
				if (sign != 0) {
					sigs[1] <- counter/nperms
				}


    cat("\nData Permutation Significance Level (for ", nperms, " permutations) = ", sigs,"\n\n\n",sep='')

	} 

npardom_output <- list(freqs=freqs, expfreqs=et, npdomfreqs=obs, 
                  expnpdomfreqs=ett, domtypes=case, kappa=kappa, 
                  z=(zkappa * sign), pk=pzkappa)

return(invisible(npardom_output))

}
