\name{bidirectional}
\alias{bidirectional}
\title{bidirectional}
\description{Tests for bidirectional dependence between pairs of lag sequential transitions.}
\usage{
bidirectional(data, labels = NULL, lag = 1, adjacent = TRUE, 
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
  \code{}Can adjacent values be coded the same? Enter "FALSE" if adjacent events 
  can never be the same. Enter "TRUE" if adjacent events can always be the same.
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
This function tests the bidirectional dependence of behaviors \emph{i} to \emph{j}, 
and \emph{j} to \emph{i}, an additive sequential pattern described by Wampold 
and Margolin (1982) and Wampold (1989, 1992). Bidirectional dependence suggests 
a reciprocal effect of behaviors. That is, behavior \emph{i} influences 
behavior \emph{j} and behavior \emph{j} influences behavior \emph{i}. For 
example, if behavior \emph{i} is a husband's positive behavior, and behavior 
\emph{j} is his wife's positive behavior, a test of bidirectional dependence 
asks whether the husband reciprocates the wife's positive behavior, \emph{and} 
the wife reciprocates the husband's positive behavior (See Margolin and Wampold, 1982). 
Bidirectional dependence is sometimes called a "circuit".
}

\value{ A list with the following elements:

  \item{freqs}{The transitional frequency matrix}

  \item{bifreqs}{The bidirectional frequencies}

  \item{expbifreqs}{The expected bidirectional frequencies}

  \item{kappas}{The bidirectional kappas}

  \item{z}{The z values for the kappas}

  \item{pk}{The p values (significance levels) for the kappas}
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
}
\author{Zakary A. Draper & Brian P. O'Connor}
\examples{
bidirectional(data_Wampold_1982, 
              labels = c('HPos','HNeu','HNeg','WPos','WNeu','WNeg'),
              permtest = TRUE, nperms = 100)
}
\keyword{ Sequential Analysis }
