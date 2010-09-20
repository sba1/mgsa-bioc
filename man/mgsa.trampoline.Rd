\name{mgsa.trampoline}
\alias{mgsa.trampoline}
\title{mgsa.trampoline}
\usage{mgsa.trampoline(o, sets, n, alpha=0.05 * (1:10), beta=0.1 * (1:8),
    p=(1:20)/length(sets), discrete=rep(TRUE, 3), alpha.breaks=NA,
    beta.breaks=NA, p.breaks=NA, steps=1e+06, restarts=1, threads=0,
    as=integer(0))}
\description{Trampoline to jump into the fast C implementation.}
\keyword{internal}
