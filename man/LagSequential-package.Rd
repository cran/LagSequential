\name{LagSequential-package}

\alias{LagSequential-package}

\docType{package}

\title{ Lag-Sequential Categorical Data Analysis }

\description{
This package provides functions for conducting lag sequential
analyses of categorical data.The functions are R versions of the programs 
provided by O'Connor (1999). The functions read a stream of codes, or a frequency
transition matrix, and produce a variety of lag sequential statistics,
including transitional frequencies, expected transitional frequencies,
transitional probabilities, z values, adjusted residuals, Yule's Q
values, likelihood ratio tests of stationarity across time and
homogeneity across groups or segments, transformed kappas for
unidirectional dependence, bidirectional dependence, parallel and
nonparallel dominance, and significance levels based on both
parametric and randomization tests.

When data is a frequency transition matrix, the code value that
preceded the first code in the sequence, and the code value that
followed the final code value, are usually unknown/unavailable.
This missing information may cause slight inaccuracies in some
of the provided statistics, most likely only at the second decimal
place. The inaccuracies will be negligible in longer data sequences.
}
\references{
 {O'Connor, B. P. (1999). Simple and flexible SAS and SPSS programs for 
 	analyzing lag-sequential categorical data. \emph{Behavior Research 
 	Methods, Instrumentation, and Computers, 31,} 718-726.}
}

