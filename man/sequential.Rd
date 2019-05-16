\name{sequential}
\alias{sequential}
\title{sequential}
\description{Computes a variety of lag sequential analysis statistics for one 
	series of codes.}
\usage{
sequential(data, labels = NULL, lag = 1, adjacent = TRUE,
           onezero = NULL, tailed = 2, permtest = FALSE, nperms = 10)
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
  a list of string variables. 
  If unspecified, codes will be labeled "Code1", "Code2", etc.
  }
  \item{lag}{
  \code{}The lag number for the analyses.
  }
  \item{adjacent}{
  \code{}Can adjacent values be coded the same? Enter "FALSE" if adjacent events 
  can never be the same. Enter "TRUE" if any adjacent events can be the 
  same. If some adjacent events can, and others cannot, be the 
  same, then enter the appropriate onezero matrix for your data using the 
  onezero argument.
  }
  \item{onezero}{
  \code{}Optional argument for specifying the one-zero matrix for the data. 
  Accepts a square matrix of ones and zeros with length ncodes. A "1" 
  indicates that the expected frequency for a given cell is to be estimated, 
  whereas a "0" indicates that the expected frequency for the cell should 
  NOT be estimated, typically because it is a structural zero (codes that 
  cannot follow one another). By default, the matrix that is created by the 
  above commands has zeros along the main diagonal, and ones everywhere 
  else, which will be appropriate for most data sets. However, if your data 
  happen to involve structural zeros that occur in cells other than the 
  cells along the main diagonal, then you must create a onezero matrix with 
  ones and zeros that is appropriate for your data.
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
Tests unidirectional dependence of states (codes). Specifically, this function
tests the hypothesis that state \emph{i} (the antecedent) follows state 
\emph{j} (the consequence) with a greater than chance probability. 
Computes a variety of statistics including two indices of effect size with 
corresponding significance tests. The larger the effect the more like the 
consequence is to follow the antecedent.
}

\value{
Displays the transitional frequency matrix, expected frequencies, transitional 
probabilities, adjusted residuals and significance levels, Yule's Q values, 
transformed Kappas (Wampold, 1989, 1992, 1995), z values for the kappas, 
and significance levels.

	Returns a list with the following elements:

  \item{freqs}{The transitional frequency matrix}

  \item{expfreqs}{The expected frequencies}

  \item{probs}{The transitional probabilities}

  \item{chi}{The overall chi-square test of the difference between the 
  	observed and expected transitional frequencies}

  \item{adjres}{The adjusted residuals}

  \item{p}{The statistical significance levels}

  \item{YulesQ}{Yule's Q values, indicating the strength of the relationships 
  	between the antecedent and the consequence transitions}

  \item{kappas}{The nonparallel dominance kappas}

  \item{z}{The z values for the kappas}

  \item{pk}{The p-values for the kappas}
}
\references{
 {O'Connor, B. P. (1999). Simple and flexible SAS and SPSS programs for analyzing 
 	lag-sequential categorical data. \emph{Behavior Research Methods, 
 		Instrumentation, and Computers, 31,} 718-726.}
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
# data is a one-column dataframe of code sequences
sequential(data_sequential, permtest = TRUE, nperms = 100)


# in this case, data is the frequency transition matrix from 
# Griffin, W. A., & Gottman, J. M. (1990). Statistical methods for analyzing family 
# interaction. In G. R. Patterson (Ed.), Family social interaction: Content and methodology
# issues in the study of aggression and depression (p. 137). Hillsdale, NJ: Erlbaum.
freqs <- t(matrix(c(
0, 0, 0, 0, 2, 2,
0,10, 5, 5,60,20,
0, 9, 2, 1, 3, 0,
0, 3, 0, 1, 5, 0,
3,54, 6, 2,24, 8,
1,24, 2, 1, 3, 12  ), 6, 6) )

sequential(freqs, adjacent = 1, 
		   labels = c('H+','Ho','H-','W+','Wo','W-'))


# Data from p 159 of Bakeman & Quera (2011), Sequential Analysis and Observational 
# Methods for the Behavioral Sciences. Cambridge University Press.
data_BQ2011 <- t(matrix(c(
2,1,4,3,3,4,3,4,2,1,4,4,5,4,1,3,4,5,3,2,2,1,4,1,2,
5,2,1,2,3,3,1,4,4,1,4,1,3,3,3,1,5,2,1,1,3,1,4,1,2,
3,3,4,5,5,2,3,3,5,2,5,4,4,2,3,1,5,5,2,2,1,3,3,3,3 )) )

sequential(data_BQ2011, labels=c('Chat','Write','Read','Ask','Attentive'),
           permtest = TRUE, nperms = 1000, tailed = 1)
}
\keyword{ Sequential Analysis }
