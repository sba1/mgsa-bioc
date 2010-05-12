\name{MgsaResults}
\alias{MgsaResults}
\title{MgsaResults}
\description{Note that representation is subject to change. Use the appropriate
accessor function to access the attributes/slots.}
\details{Here are some internal details of the current implementation.

Attributes alphaPost, betaPost, and pPost are data frames, of which the first column
represents a realization value, and subsequent columns the posterior probabilities of
each MCMC run for that realization.

Attribute setsResults is a dataFrame in which the first columns represent the marginal
posteriors of each set. The last two hold, for each set, the counts.

The current representation is tightly coupled to an MCMC solver, which may
or may not replaced in the future. Therefore, everything is subject to change.
You've been warned!}
\seealso{\code{\link{populationSize}}
\code{\link{studySetSizeInPopulation}}
\code{\link{alphaPost}}
\code{\link{betaPost}}
\code{\link{pPost}}}
