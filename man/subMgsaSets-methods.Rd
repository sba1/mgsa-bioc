\name{subMgsaSets-methods}
\alias{subMgsaSets-methods}
\alias{subMgsaSets}
\title{Subset of an MgsaSets}
\description{Returns a subset of an \code{\linkS4class{MgsaSets}} in which only the required
items are kept. Empty sets are removed.}
\note{TODO: do we need that method? Does it work on items or item indexes?}
\value{\code{subMgsaSets-methods}: an \code{\linkS4class{MgsaSets}}.

}
\arguments{\item{sets}{an \code{\linkS4class{MgsaSets}}.}
\item{items}{\code{numeric}. The items to restrict on.}
}
\alias{subMgsaSets,MgsaSets,numeric-method}
\alias{subMgsaSets}
\usage{\S4method{subMgsaSets}{MgsaSets,numeric}(sets, items)
}

