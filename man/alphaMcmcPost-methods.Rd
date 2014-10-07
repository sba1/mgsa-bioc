\name{alphaMcmcPost}
\alias{alphaMcmcPost}
\alias{alphaMcmcPost-methods}
\alias{alphaMcmcPost,MgsaMcmcResults-method}
\title{posterior estimates of the parameter alpha for each MCMC run}
\usage{
  alphaMcmcPost(x)

  \S4method{alphaMcmcPost}{MgsaMcmcResults}(x)
}
\arguments{
  \item{x}{a \code{\linkS4class{MgsaMcmcResults}}.}
}
\value{
  \code{matrix}: Posterior estimates of the parameter alpha
  for each MCMC run.
}
\description{
  Posterior estimates of the parameter alpha for each MCMC
  run.
}

