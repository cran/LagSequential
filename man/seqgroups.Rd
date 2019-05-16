\name{seqgroups}
\alias{seqgroups}
\title{seqgroups}
\description{Computes a variety of sequential analysis statistics 
	for data that are in segments (e.g, for multiple dyads or groups).}
\usage{
seqgroups(alldata, labels = NULL, lag = 1, adjacent = TRUE,
          onezero = NULL, tailed = 2, test = "homogeneity", 
          output = "all")
}
\arguments{
  \item{alldata}{
  \code{}A one-column dataframe, or a vector of code sequences, which
  can be numeric or strings. Missing values are not permitted.\cr
  
  If alldata is numeric, then the integers must range from "1" to the total number of 
  possible code values (which is not the total number of codes in a sequence), 
  and a number greater than 999 must be used in alldata to
  separate the codes sequences for different groups/dyads.
  See "data_seqgroups" for an example.\cr
  
  If alldata consists of strings/characters, then the
  word "segment" must be used in alldata to separate the code sequences
  for different groups/dyads. See "data_seqgroups" for an example.
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
  can never be the same. Enter "TRUE" if adjacent events can always be the same. 
  Enter "TRUE" if some adjacent events can, and others cannot, be the same; 
  then enter the appropriate onezero matrix for your data.
  }
  \item{onezero}{
  \code{}Optional argument for specifying the one-zero matrix for the data.
  Useful when some adjacent events can, and others cannot, be the same.
  Accepts a square matrix of ones and zeros with length ncodes. A "1" indicates 
  that the expected frequency for a given cell is to be estimated, whereas a 
  "0" indicates that the expected frequency for the cell should NOT be estimated, 
  typically because it is a structural zero (codes that cannot follow one another). 
  By default, the matrix that is created by the above commands has zeros along the 
  main diagonal, and ones everywhere else, which will be appropriate for most data 
  sets. However, if your data happen to involve structural zeros that occur in cells 
  other than the cells along the main diagonal, then you must create a onezero 
  matrix with ones and zeros that is appropriate for your data.
  }
  \item{tailed}{
  \code{}Specify whether significance tests are one-tailed or two-tailed. 
  Options are "1" or "2".
  }
  \item{test}{
  \code{}Specify whether to run tests for homogeneity of homogeneity or 
  stationarity. Homogeneity should be tested when groups in the data are 
  actually different groups, whereas stationarity should be tested when 
  groups in the data are segments of a single stream of observations. 
  Options are "homogeneity" or "stationarity".
  }
  \item{output}{
  \code{}Specify the desired output. Options are "pooled" for pooled data only, 
  or "all" for all data sets.
  }
}
\details{
Computes a variety of sequential analysis statistics for data that are 
in segments (e.g, for multiple dyads or groups. This is the same as the 
"sequential" function provided in this package, but allows for the data 
to be segmented. Sequential statistics are calculated for each segment, 
as well as for the data pooled across all segments.
}

\value{
For each of the groups or segments and for the pooled data, displays the 
transitional frequency matrix, expected frequencies, transitional probabilities, 
adjusted residuals and significance levels, Yule's Q values, transformed 
Kappas (Wampold, 1989, 1992, 1995), z values for the kappas, and significance levels.

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

  \item{output}{The requested output data}
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
seqgroups(data_seqgroups_strings)
}
\keyword{ Sequential Analysis }
