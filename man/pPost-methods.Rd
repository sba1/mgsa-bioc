\name{pPost}
\alias{pPost}
\alias{pPost-methods}
\alias{pPost,MgsaResults-method}
\title{Posterior for beta}
\usage{
  pPost(x)

  \S4method{pPost}{MgsaResults}(x)
}
\arguments{
  \item{x}{a \code{\linkS4class{MgsaResults}}.}
}
\value{
  \code{data.frame}: realization values, posterior estimate
  and standard error for the parameter p.
}
\description{
  Realization values, posterior estimate and standard error
  for the parameter p.
}

