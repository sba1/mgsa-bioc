\name{mgsa}
\alias{mgsa}
\alias{mgsa,character,list-method}
\alias{mgsa,character,MgsaSets-method}
\alias{mgsa,integer,list-method}
\alias{mgsa,logical,list-method}
\alias{mgsa-methods}
\alias{mgsa,numeric,list-method}
\title{Performs an MGSA analysis}
\usage{
  mgsa( o, sets, population=NULL, alpha=seq(0.01,0.3,
    length.out=10), beta=seq(0.1,0.8, length.out=10),
    p=seq( min(0.1, 1/length(sets)), min(0.3,
    20/length(sets)), length.out=10), steps=1e6,
    restarts=5, threads=0 )

  \S4method{mgsa}{integer,list}(o, sets, population = NULL,
    alpha = seq(0.01, 0.3, length.out = 10),
    beta = seq(0.1, 0.8, length.out = 10),
    p = seq(min(0.1, 1/length(sets)), min(0.3, 20/length(sets)), length.out = 10),
    steps = 1e+06, restarts = 5, threads = 0)

  \S4method{mgsa}{numeric,list}(o, sets, population = NULL,
    alpha = seq(0.01, 0.3, length.out = 10),
    beta = seq(0.1, 0.8, length.out = 10),
    p = seq(min(0.1, 1/length(sets)), min(0.3, 20/length(sets)), length.out = 10),
    steps = 1e+06, restarts = 5, threads = 0)

  \S4method{mgsa}{character,list}(o, sets,
    population = NULL,
    alpha = seq(0.01, 0.3, length.out = 10),
    beta = seq(0.1, 0.8, length.out = 10),
    p = seq(min(0.1, 1/length(sets)), min(0.3, 20/length(sets)), length.out = 10),
    steps = 1e+06, restarts = 5, threads = 0)

  \S4method{mgsa}{logical,list}(o, sets, population = NULL,
    alpha = seq(0.01, 0.3, length.out = 10),
    beta = seq(0.1, 0.8, length.out = 10),
    p = seq(min(0.1, 1/length(sets)), min(0.3, 20/length(sets)), length.out = 10),
    steps = 1e+06, restarts = 5, threads = 0)

  \S4method{mgsa}{character,MgsaSets}(o, sets,
    population = NULL,
    alpha = seq(0.01, 0.3, length.out = 10),
    beta = seq(0.1, 0.8, length.out = 10),
    p = seq(min(0.1, 1/length(sets)), min(0.3, 20/length(sets)), length.out = 10),
    steps = 1e+06, restarts = 5, threads = 0)
}
\arguments{
  \item{o}{The observations. It can be a \code{numeric},
  \code{integer}, \code{character} or \code{logical}. See
  details.}

  \item{sets}{The sets. It can be an
  \code{\linkS4class{MgsaSets}} or a \code{list}. In this
  case, each list entry is a vector of type \code{numeric},
  \code{integer}, \code{character}. See details.}

  \item{population}{The total population. Optional. A
  \code{numeric}, \code{integer} or \code{character}
  vector. Default to \code{NULL}. See details.}

  \item{alpha}{Grid of values for the parameter alpha.
  Values represent probabilities of false-positive events
  and hence must be in [0,1]. \code{numeric}.}

  \item{beta}{Grid of values for the parameter beta. Values
  represent probabilities of false-negative events and
  hence must be in [0,1]. \code{numeric}.}

  \item{p}{Grid of values for the parameter p. Values
  represent probabilities of term activity and therefore
  must be in [0,1]. \code{numeric}.}

  \item{steps}{The number of steps of each run of the MCMC
  sampler. \code{integer} of length 1. A recommended value
  is 1e6 or greater.}

  \item{restarts}{The number of different runs of the MCMC
  sampler. \code{integer} of length 1. Must be greater or
  equal to 1. A recommended value is 5 or greater.}

  \item{threads}{The number of threads that should be used
  for concurrent restarts. A value of 0 means to use all
  available cores. Default to 0.}
}
\value{
  An \code{\link{MgsaMcmcResults}} object.
}
\description{
  Estimate marginal posterior of the MGSA problem with an
  MCMC sampling algorithm.
}
\details{
  The function can handle items (such as genes) encoded as
  \code{character} or \code{integer}. For convenience
  \code{numeric} items can also be provided but these
  values should essentially be integers. The type of items
  in the observations \code{o}, the \code{sets} and in the
  optional \code{population} should be consistent. In the
  case of \code{character} items, \code{o} and
  \code{population} should be of type \code{character} and
  \code{sets} can either be an
  \code{\linkS4class{MgsaSets}} or a \code{list} of
  \code{character} vectors. In the case of \code{integer}
  items, \code{o} should be of type \code{integer},
  \code{numeric} (but essentially with integer values), or
  \code{logical} and entries in \code{sets} as well as the
  \code{population} should be \code{integer}. When \code{o}
  is \code{logical}, it is first coerced to integer with a
  call on \code{\link{which}}. Observations outside the
  \code{population} are not taken into account. If
  \code{population} is \code{NULL}, it is defined as the
  union of all sets.

  The default grid value for p is such that between 1 and
  20 sets are active in expectation. The lower limit is
  constrained to be lower than 0\.1 and the upper limit
  lower than 0\.3 independently of the total number of sets
  to make sure that complex solutions are penalized.
  Marginal posteriors of activity of each set are estimated
  using an MCMC sampler as described in Bauer et al., 2010.
  Because convergence of an MCM sampler is difficult to
  assess, it is recommended to run it several times (using
  \code{restarts}). If variations between runs are too
  large (see \code{\link{MgsaResults}}), the number of
  steps (\code{steps}) of each MCMC run should be
  increased.
}
\examples{
## observing items A and B, with sets {A,B,C} and {B,C,D}
mgsa(c("A", "B"), list(set1 = LETTERS[1:3], set2 = LETTERS[2:4]))

## same case with integer representation of the items and logical observation
mgsa(c(TRUE,TRUE,FALSE,FALSE), list(set1 = 1:3, set2 = 2:4))

## a small example with gene ontology sets and plot
data(example)
fit = mgsa(example_o, example_go)
## Not run:
plot(fit)
## End(Not run)
}
\references{
  Bauer S., Gagneur J. and Robinson P. GOing Bayesian:
  model-based gene set analysis of genome-scale data.
  Nucleic Acids Research (2010)
  \url{http://nar.oxfordjournals.org/content/38/11/3523.full}
}
\seealso{
  \code{\link{MgsaResults}}, \code{\link{MgsaMcmcResults}}
}

