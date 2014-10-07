\name{itemIndices}
\alias{itemIndices}
\alias{itemIndices-methods}
\alias{itemIndices,MgsaSets,character-method}
\alias{itemIndices,MgsaSets,numeric-method}
\title{Item indices of a MgsaSets}
\usage{
  itemIndices(sets, items)

  \S4method{itemIndices}{MgsaSets,character}(sets, items)

  \S4method{itemIndices}{MgsaSets,numeric}(sets, items)
}
\arguments{
  \item{sets}{an instance of class
  \code{\linkS4class{MgsaSets}}.}

  \item{items}{\code{character} or \code{numeric} the items
  of interest.}
}
\value{
  a \code{integer}: the item indices.
}
\description{
  Returns the indices corresponding to the items
}

