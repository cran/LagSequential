\name{nonparadom}
\alias{nonparadom}
\title{nonparadom}
\description{Tests for nonparallel dominance, a form of asymmetry in 
	predictability, between i to j and k to L (Wampold, 1984, 1989, 1992, 1995).}
\usage{
nonparadom(data, i, j, k, L, labels = NULL, lag = 1, adjacent = TRUE, 
           tailed = 1, permtest = FALSE, nperms = 10)
}
\arguments{
  \item{data}{
  \code{}A one-column dataframe, or a vector of code sequences, or a square
  frequency transition matrix. If data is not a frequency transition matrix,
  then data must be either (a) a series of string (non-numeric) code values,
  or (b) a series of integer codes with values ranging from "1" to what 
  ever value the user specifies in the "ncodes" argument. There should be no
  code values with zero frequencies. Missing values are not permitted.
}
  \item{i}{
  \code{}Code value for i.
  }
  \item{j}{
  \code{}Code value for j.
  }
  \item{k}{
  \code{}Code value for k.
  }
  \item{L}{
  \code{}Code value for L.
  }
  \item{labels}{
  \code{}Optional argument for providing labels to the code values. Accepts a 
  list of string variables. If 
  unspecified, codes will be labeled "Code1", "Code2", etc.
  }
  \item{lag}{
  \code{}The lag number for the analyses.
  }
  \item{adjacent}{
  \code{}Can adjacent values be coded the same? Options are "TRUE" for yes or "FALSE" for no.
  }
  \item{tailed}{
  \code{}Specify whether significance tests are one-tailed or two-tailed. Options are "1" or "2".
  }
  \item{permtest}{
  \code{}Do you want to run permutation tests of significance? 
  Options are "FALSE" for no, or "TRUE" for yes. Warning: these computations can be time consuming.
  }
  \item{nperms}{
  \code{}The number of permutations per block.
  }
}
\details{
Tests for nonparallel dominance or asymmetry in predictability, which is the 
difference in predictability between i to j and k to L, as described by 
Wampold (1984, 1989, 1992, 1995).
Parallel dominance (another function in this package) is the difference in 
predictability between i to j and j to i.
In parallel dominance the i and j values across the two pairs of codes
are the same. In nonparallel dominance, the i and j values across the
two pairs of codes may vary, i.e., they do not have to be the same.
}

\value{Displays the transitional frequency matrix, expected frequencies, 
	expected and observed nonparallel dominance frequencies, kappas, 
	the z values for the kappas, and the significance levels.
	
	Returns a list with the following elements:

  \item{freqs}{The transitional frequency matrix}

  \item{expfreqs}{The expected frequencies}

  \item{npdomfreqs}{The nonparallel dominance frequencies}

  \item{expnpdomfreqs}{The expected nonparallel dominance frequencies}

  \item{domtypes}{There are 4 sequential dominance case types 
  	described by Wampold (1989). These cases describe the direction 
  	of the effect for \emph{i} on \emph{j} and \emph{j} on \emph{i}. 
  	The four cases are: (1) \emph{i} increases \emph{j}, and \emph{j} 
  	increases \emph{i}, (2) \emph{i} decreases \emph{j}, and \emph{j} 
  	decreases \emph{i}, (3) \emph{i} increases \emph{j}, and \emph{j} 
  	decreases \emph{i}, and (4) \emph{i} decreases \emph{j}, and \emph{j} 
  	increases \emph{i}. Each cell of this matrix indicates the case that 
  	applies to the transition indicated by the cell.}

  \item{kappas}{The nonparallel dominance kappas}

  \item{z}{The z values for the kappas}

  \item{pk}{The p-values for the kappas}
}

\references{
 {O'Connor, B. P. (1999). Simple and flexible SAS and SPSS programs for analyzing 
 	lag-sequential categorical data. \emph{Behavior Research Methods, 
 		Instrumentation, and Computers, 31,} 718-726.}
\cr\cr {Wampold, B. E., & Margolin, G. (1982). Nonparametric strategies to test 
	the independence of behavioral states in sequential data. \emph{Psychological 
		Bulletin, 92,} 755-765.}
\cr\cr {Wampold, B. E. (1984). Tests of dominance in sequential categorical data. 
	\emph{Psychological Bulletin, 96,} 424-429.}
\cr\cr {Wampold, B. E. (1989). Kappa as a measure of pattern in sequential data. 
	\emph{Quality & Quantity, 23,} 171-187.}
\cr\cr {Wampold, B. E. (1992). The intensive examination of social interactions. 
	In T. Kratochwill & J. Levin (Eds.), \emph{Single-case research design and 
		analysis: New directions for psychology and education} (pp. 93-131). 
		Hillsdale, NJ: Erlbaum.}
\cr\cr {Wampold, B. E. (1995). Analysis of behavior sequences in psychotherapy. 
	In J. Siegfried (Ed.), \emph{Therapeutic and everyday discourse as behavior 
		change: Towards a micro-analysis in psychotherapy process research} 
		(pp. 189-214). Norwood, NJ: Ablex.}
}
\author{Zakary A. Draper & Brian P. O'Connor}
\examples{
nonparadom(data_Wampold_1984, i = 6, j = 1, k = 3, L = 4,
           labels = c('HPos','HNeu','HNeg','WPos','WNeu','WNeg'), 
           permtest = TRUE, nperms = 1000)
}
\keyword{ Sequential Analysis }
