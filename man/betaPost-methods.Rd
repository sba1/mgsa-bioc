\name{betaPost}
\alias{betaPost}
\alias{betaPost-methods}
\alias{betaPost,MgsaResults-method}
\title{Posterior for beta}
\usage{
  betaPost(x)

  \S4method{betaPost}{MgsaResults}(x)
}
\arguments{
  \item{x}{a \code{\linkS4class{MgsaResults}}.}
}
\value{
  \code{data.frame}: realization values, posterior estimate
  and standard error for the parameter beta.
}
\description{
  Realization values, posterior estimate and standard error
  for the parameter beta.
}

