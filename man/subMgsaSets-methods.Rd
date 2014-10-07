\name{subMgsaSets}
\alias{subMgsaSets}
\alias{subMgsaSets-methods}
\alias{subMgsaSets,MgsaSets,character-method}
\title{Subset of an MgsaSets}
\usage{
  subMgsaSets(sets, items)

  \S4method{subMgsaSets}{MgsaSets,character}(sets, items)
}
\arguments{
  \item{sets}{an \code{\linkS4class{MgsaSets}}.}

  \item{items}{\code{character}. The items to restrict on.}
}
\value{
  an \code{\linkS4class{MgsaSets}}.
}
\description{
  Returns a subset of an \code{\linkS4class{MgsaSets}} that
  contains only the specified items. Empty sets are
  removed.
}

