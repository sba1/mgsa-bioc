\name{setsResults}
\alias{setsResults}
\alias{setsResults-methods}
\alias{setsResults,MgsaResults-method}
\title{Posterior for each set}
\usage{
  setsResults(x)

  \S4method{setsResults}{MgsaResults}(x)
}
\arguments{
  \item{x}{a \code{\linkS4class{MgsaResults}}.}
}
\value{
  \code{data.frame}: For each set, number of items of the
  set in the population, number of items of the set in the
  study set, estimate of its marginal posterior probability
  and standard error.
}
\description{
  Number of items of the set in the population, the number
  of items of the set in the study set, the estimate of its
  marginal posterior probability and its standard error.
}

