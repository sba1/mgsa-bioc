\name{alphaPost}
\alias{alphaPost}
\alias{alphaPost-methods}
\alias{alphaPost,MgsaResults-method}
\title{Posterior for alpha}
\usage{
  alphaPost(x)

  \S4method{alphaPost}{MgsaResults}(x)
}
\arguments{
  \item{x}{a \code{\linkS4class{MgsaResults}}.}
}
\value{
  \code{data.frame}: realization values, posterior estimate
  and standard error for the parameter alpha.
}
\description{
  Realization values, posterior estimate and standard error
  for the parameter alpha.
}

