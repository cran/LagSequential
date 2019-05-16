\name{twocells}
\alias{twocells}
\title{twocells}
\description{Simultaneously tests the unidirectional dependence of 
	\emph{i} to \emph{j}, and the unidirectional dependence of \emph{k} 
	to \emph{L}, an additive pattern described by Wampold and Margolin 
	(1982) and Wampold (1989, 1992).}
\usage{
twocells(data, i, j, k, L, labels = NULL, lag = 1,
         adjacent = TRUE, tailed = 1, permtest = FALSE, nperms = 10)
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
  \code{}Code value for \emph{i}.
  }
  \item{j}{
  \code{}Code value for \emph{j}.
  }
  \item{k}{
  \code{}Code value for \emph{k}.
  }
  \item{L}{
  \code{}Code value for \emph{L}.
  }
  \item{labels}{
  \code{}Optional argument for providing labels to the code values. Accepts 
  a list of string variables. 
  If unspecified, codes will be labeled "Code1", "Code2", etc.
  }
  \item{lag}{
  \code{}The lag number for the analyses.
  }
  \item{adjacent}{
  \code{}Can adjacent values be coded the same? Options are "TRUE" for yes, 
  and "FALSE" for no.
  }
  \item{tailed}{
  \code{}Specify whether significance tests are one-tailed or two-tailed. 
  Options are "1" or "2".
  }
  \item{permtest}{
  \code{}Do you want to run permutation tests of significance? Options 
  are "FALSE" for no, or "TRUE" for yes. Warning: these computations can be time consuming.
  }
  \item{nperms}{
  \code{}The number of permutations per block.
  }
}
\details{
This function simultaneously tests the unidirectional dependence of \emph{i} 
to \emph{j} and the unidirectional dependence of \emph{k} to \emph{L}. The 
user specifies the code values used for \emph{i}, \emph{j}, \emph{k}, and 
\emph{L} in the analyses. For example, Wampold and Margolin (1982) 
described a situation wherein a spouse responds to negative behaviors 
with something other than a negative behavior.
}

\value{
Displays the transitional frequency matrix, observed and expected values for 
the two cell test, kappa, the z value for kappa, and the significance level.

	Returns a list with the following elements:

  \item{freqs}{The transitional frequency matrix}

  \item{expfreqs}{The expected frequencies}

  \item{twocellfreq}{The observed number of transitions from \emph{i} 
  	to \emph{j} and from \emph{k} to \emph{L}.}

  \item{kappa}{The twocells kappa}

  \item{z}{The z value for the kappa}

  \item{pk}{The p-value for the kappa}
}
\references{
 {O'Connor, B. P. (1999). Simple and flexible SAS and SPSS programs for analyzing 
 	lag-sequential categorical data. \emph{Behavior Research Methods, 
 		Instrumentation, and Computers, 31,} 718-726.}
\cr\cr {Wampold, B. E., & Margolin, G. (1982). Nonparametric strategies to test 
	the independence of behavioral states in sequential data. \emph{Psychological 
		Bulletin, 92,} 755-765.}
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
twocells(data_Wampold_1982, i = 6, j = 1, k = 3, L = 4,
         labels = c('HPos','HNeu','HNeg','WPos','WNeu','WNeg'),
         permtest = TRUE, nperms = 100)
}
\keyword{ Sequential Analysis }
