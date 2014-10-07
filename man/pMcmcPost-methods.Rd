\name{pMcmcPost}
\alias{pMcmcPost}
\alias{pMcmcPost-methods}
\alias{pMcmcPost,MgsaMcmcResults-method}
\title{posterior estimates of the parameter p for each MCMC run}
\usage{
  pMcmcPost(x)

  \S4method{pMcmcPost}{MgsaMcmcResults}(x)
}
\arguments{
  \item{x}{a \code{\linkS4class{MgsaMcmcResults}}.}
}
\value{
  \code{matrix}: Posterior estimates of the parameter p for
  each MCMC run.
}
\description{
  Posterior estimates of the parameter p for each MCMC run.
}

