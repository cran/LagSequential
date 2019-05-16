\name{paradom}
\alias{paradom}
\title{paradom}
\description{Tests for parallel dominance in lag sequential data.}
\usage{
paradom(data, labels = NULL, lag = 1, adjacent = TRUE,
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
  \item{labels}{
  \code{}Optional argument for providing labels to the code values. Accepts 
  a list of string variables. If 
  unspecified, codes will be labeled "Code1", "Code2", etc.
  }
  \item{lag}{
  \code{}The lag number for the analyses.
  }
  \item{adjacent}{
  \code{}Can adjacent values be coded the same? Options are "TRUE" for yes or "FALSE" for no.
  }
  \item{tailed}{
  \code{}Specify whether significance tests are one-tailed or two-tailed. 
  Options are "1" or "2".
  }
  \item{permtest}{
  \code{}Do you want to run permutation tests of significance? Options are 
  "FALSE" for no, or "TRUE" for yes. Warning: these computations can be time consuming.
  }
  \item{nperms}{
  \code{}The number of permutations per block.
  }
}
\details{
Tests for parallel dominance or asymmetry in predictability, which is 
the difference in predictability between i to j and j to i (e.g., whether 
B's behavior is more predictable from A's behavior than vice versa), 
as described by Wampold (1984, 1989, 1992, 1995).
}
\value{
Displays the transitional frequency matrix and matrices of expected frequencies, 
expected and observed parallel dominance frequencies, parallel dominance kappas, 
z values for the kappas, and significance levels. There are four possible cases, 
or kinds, of parallel dominance (see Wampold 1989, 1992, 1995), and the function 
returns a matrix indicating the kind of case for each cell in the transitional 
frequency matrix.

	Returns a list with the following elements:

  \item{freqs}{The transitional frequency matrix}

  \item{expfreqs}{The expected frequencies}

  \item{domfreqs}{The parallel dominance frequencies}

  \item{expdomfreqs}{The expected parallel dominance frequencies}

  \item{domtypes}{There are 4 sequential dominance case types described 
  	by Wampold (1989). These cases describe the direction of the effect for 
  	\emph{i} on \emph{j} and \emph{j} on \emph{i}. The four cases are: (1) 
  	\emph{i} increases \emph{j}, and \emph{j} increases \emph{i}, (2) \emph{i} 
  	decreases \emph{j}, and \emph{j} decreases \emph{i}, (3) \emph{i} increases 
  	\emph{j}, and \emph{j} decreases \emph{i}, and (4) \emph{i} decreases \emph{j}, 
  	and \emph{j} increases \emph{i}. Each cell of this matrix indicates the case 
  	that applies to the transition indicated by the cell.}

  \item{kappas}{The parallel dominance kappas}

  \item{z}{The z values for the kappas}

  \item{pk}{The p-values for the kappas}
}
\references{
 {O'Connor, B. P. (1999). Simple and flexible SAS and SPSS programs for analyzing 
 	lag-sequential categorical data. \emph{Behavior Research Methods, 
 		Instrumentation, and Computers, 31,} 718-726.}
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
paradom(data_Wampold_1984, 
        labels = c('HPos','HNeu','HNeg','WPos','WNeu','WNeg'), 
        permtest = TRUE, nperms = 1000)
}
\keyword{ Sequential Analysis }
