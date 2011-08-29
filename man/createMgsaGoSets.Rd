\name{createMgsaGoSets}
\alias{createMgsaGoSets}
\title{This functions takes a 1:1 mapping of go...}
\usage{createMgsaGoSets(go.ids, items)
}
\description{This functions takes a 1:1 mapping of go.ids to items and returns
a full MgsaGOSets instance. The structure of GO is gathered from GO.db. It
is sufficient to specify just the directly asserted mapping (or annotation), 
i.e., the most specific ones. The true path rule is taken account, that is, if an
item is annotated to a term then it will be also annotated to more general
terms (some people prefer to say that just the transitive closure is calculated).}
\arguments{\item{go.ids}{a character vector of GO ids (GO:00001234)}
\item{items}{a vector of identifiers that are annotated to the term
in the corresponding position of the go.ids vector.}
}

