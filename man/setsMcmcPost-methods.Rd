\name{setsMcmcPost}
\alias{setsMcmcPost}
\alias{setsMcmcPost-methods}
\alias{setsMcmcPost,MgsaMcmcResults-method}
\title{posterior estimates of the the set marginal probabilities  for each MCMC run}
\usage{
  setsMcmcPost(x)

  \S4method{setsMcmcPost}{MgsaMcmcResults}(x)
}
\arguments{
  \item{x}{a \code{\linkS4class{MgsaMcmcResults}}.}
}
\value{
  \code{matrix}: Posterior estimates of the set marginal
  probabilities for each MCMC run.
}
\description{
  Posterior estimates of the set marginal probabilities for
  each MCMC run.
}

