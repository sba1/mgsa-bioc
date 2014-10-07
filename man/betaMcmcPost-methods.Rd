\name{betaMcmcPost}
\alias{betaMcmcPost}
\alias{betaMcmcPost-methods}
\alias{betaMcmcPost,MgsaMcmcResults-method}
\title{posterior estimates of the parameter beta for each MCMC run}
\usage{
  betaMcmcPost(x)

  \S4method{betaMcmcPost}{MgsaMcmcResults}(x)
}
\arguments{
  \item{x}{a \code{\linkS4class{MgsaMcmcResults}}.}
}
\value{
  \code{matrix}: Posterior estimates of the parameter beta
  for each MCMC run.
}
\description{
  Posterior estimates of the parameter beta for each MCMC
  run.
}

