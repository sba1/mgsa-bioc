## the class of MgsaResults: Class of the returned value of mgsa()
## Useful for plotting and show purposes

######## Class definitions
#### MgsaResults

#'
#' The results of an MGSA analysis.
#' 
#' @slot populationSize The number of items in the population. 
#' @slot studySetSizeInPopulation The number of items both in the study set and in the population.
#' @slot alphaPost with columns \code{value}, \code{estimate} and \code{std.error}.
#' @slot betaPost with columns \code{value}, \code{estimate} and \code{std.error}.
#' @slot pPost with columns \code{value}, \code{estimate} and \code{std.error}.
#' @slot setsResults with columns \code{inPopulation}, \code{inStudySet}, \code{estimate} and \code{std.error}. 
#'   
#' The columns of the slots \code{alphaPost}, \code{betaPost}, and \code{pPost} contains a realization value, its posterior estimate and standard error for the parameters alpha, beta and p respectively.
#' 
#' The columns of the slot \code{setsResults} contains the number of items of the set in the population, the number of items of the set in the study set, the estimate of its marginal posterior probability and its standard error.
#' The \code{\link{rownames}} are the names of the sets if available.
#' 
#' Because an \code{MgsaResults} is the outcome of an MGSA analysis (see \code{\link{mgsa}}), accessors but no replacement methods exist for each slot.
#' 
#' @title Results of an MGSA analysis
#' @seealso \code{\link{mgsa}}
#' @exportClass MgsaResults

setClass(
		"MgsaResults",
		representation = representation(
				populationSize = "numeric",
				studySetSizeInPopulation="numeric",
				alphaPost = "data.frame", betaPost ="data.frame", pPost ="data.frame",
				setsResults = "data.frame"
		)
)

### populationSize
#' The size of the population on which the analysis was run.
#' 
#' @title Size of the population of a MgsaResults
#' @param x a \code{\linkS4class{MgsaResults}}.
#' @return \code{integer}: the size of the population.
#' @rdname populationSize-methods
#' @exportMethod populationSize
setGeneric( "populationSize", function(x) standardGeneric( "populationSize" ) )


#' @rdname populationSize-methods
setMethod(
		"populationSize",
		signature=c( "MgsaResults" ),
		function( x ) x@populationSize
)

### studySetSizeInPopulation
#' The size of the study set on which the analysis was run.
#' 
#' @title Size of the study set of a MgsaResults
#' @param x a \code{\linkS4class{MgsaResults}}.
#' @return \code{integer}: the size of the study set.
#' @rdname studySetSizeInPopulation-methods
#' @exportMethod studySetSizeInPopulation
setGeneric( "studySetSizeInPopulation", function(x) standardGeneric( "studySetSizeInPopulation" ) )

#' @rdname studySetSizeInPopulation-methods
setMethod(
		"studySetSizeInPopulation",
		signature=c( "MgsaResults" ),
		function( x ) x@studySetSizeInPopulation
)


### alphaPost
#' Realization values, posterior estimate and standard error for the parameter alpha.
#' 
#' @title Posterior for alpha
#' @param x a \code{\linkS4class{MgsaResults}}.
#' @return \code{data.frame}: realization values, posterior estimate and standard error for the parameter alpha.
#' @rdname alphaPost-methods
#' @exportMethod alphaPost
setGeneric( "alphaPost", function(x) standardGeneric( "alphaPost" ) )

#' @rdname alphaPost-methods
setMethod(
		"alphaPost",
		signature=c( "MgsaResults" ),
		function( x ) x@alphaPost
)

### betaPost
#' Realization values, posterior estimate and standard error for the parameter beta.
#' 
#' @title Posterior for beta
#' @param x a \code{\linkS4class{MgsaResults}}.
#' @return \code{data.frame}: realization values, posterior estimate and standard error for the parameter beta.
#' @rdname betaPost-methods
#' @exportMethod betaPost
setGeneric( "betaPost", function(x) standardGeneric( "betaPost" ) )

#' @rdname betaPost-methods
setMethod(
		"betaPost",
		signature=c( "MgsaResults" ),
		function( x ) x@betaPost
)

### pPost
#' Realization values, posterior estimate and standard error for the parameter p.
#' 
#' @title Posterior for beta
#' @param x a \code{\linkS4class{MgsaResults}}.
#' @return \code{data.frame}: realization values, posterior estimate and standard error for the parameter p.
#' @rdname pPost-methods
#' @exportMethod pPost
setGeneric( "pPost", function(x) standardGeneric( "pPost" ) )

#' @rdname pPost-methods
setMethod(
		"pPost",
		signature=c( "MgsaResults" ),
		function( x ) x@pPost
)

### setsResults
#' Number of items of the set in the population, the number of items of the set in the study set, the estimate of its marginal posterior probability and its standard error.
#' 
#' @title Posterior for each set
#' @param x a \code{\linkS4class{MgsaResults}}.
#' @return \code{data.frame}: For each set, number of items of the set in the population, number of items of the set in the study set, estimate of its marginal posterior probability and standard error.
#' @rdname setsResults-methods
#' @exportMethod setsResults
setGeneric( "setsResults", function(x) standardGeneric( "setsResults" ) )

#' @rdname setsResults-methods
setMethod(
		"setsResults",
		signature=c( "MgsaResults" ),
		function( x ) x@setsResults
)


#'
#' Instances of this class are used to hold the additional information 
#' that was provided by running (possibly multiple times) an MCMC algorithm.
#' 
#' @slot steps how many steps per MCMC run
#' @slot restarts how many MCMC runs
#' @slot alphaMcmcPost posterior estimates for each MCMC run of the parameter alpha
#' @slot betaMcmcPost posterior estimates for each MCMC run of the parameter beta
#' @slot pMcmcPost posterior estimates for each MCMC run of the parameter p
#' @slot setsMcmcPost posterior estimates for each MCMC run of the sets marginal posterior probabilities
#' 
#' The columns of the matrices \code{alphaMcmcPost}, \code{betaMcmcPost}, \code{pMcmcPost} and setsMcmcPost stores the posterior estimates for each individual MCMC run.
#' The row order matches the one of the slot \code{alphaPost}, \code{betaPost}, \code{pPots}, and \code{setsResults} respectively.
#'  
#' Accessor methods exist for each slot.
#' 
#' @seealso \code{\link{mgsa}}
#' @exportClass MgsaMcmcResults
setClass(
		"MgsaMcmcResults",
		contains = c("MgsaResults"),
		representation = representation(
				steps = "numeric",
				restarts = "numeric",
				alphaMcmcPost = "matrix", betaMcmcPost = "matrix", pMcmcPost ="matrix", setsMcmcPost = "matrix"
		)
)

#### Number of steps
#' how many steps per MCMC run.
#' 
#' @title How many steps per MCMC run
#' @param x a \code{\linkS4class{MgsaMcmcResults}}.
#' @return \code{integer}: how many steps per MCMC run.
#' @rdname steps-methods
#' @exportMethod steps
setGeneric( "steps", function(x) standardGeneric( "steps" ) )

#' @rdname steps-methods
setMethod(
		"steps",
		signature=c( "MgsaMcmcResults" ),
		function( x ) x@steps
)

#### Number of restarts
#' how many MCMC runs.
#' 
#' @title How many MCMC runs
#' @param x a \code{\linkS4class{MgsaMcmcResults}}.
#' @return \code{integer}: how many MCMC runs.
#' @rdname restarts-methods
#' @exportMethod restarts
setGeneric( "restarts", function(x) standardGeneric( "restarts" ) )

#' @rdname restarts-methods
setMethod(
		"restarts",
		signature=c( "MgsaMcmcResults" ),
		function( x ) x@restarts
)

#### alphaMcmcPost 
#' Posterior estimates of the parameter alpha for each MCMC run.
#' 
#' @title posterior estimates of the parameter alpha for each MCMC run
#' @param x a \code{\linkS4class{MgsaMcmcResults}}.
#' @return \code{matrix}: Posterior estimates of the parameter alpha for each MCMC run.
#' @rdname alphaMcmcPost-methods
#' @exportMethod alphaMcmcPost
setGeneric( "alphaMcmcPost", function(x) standardGeneric( "alphaMcmcPost" ) )

#' @rdname alphaMcmcPost-methods
setMethod(
		"alphaMcmcPost",
		signature=c( "MgsaMcmcResults" ),
		function( x ) x@alphaMcmcPost
)


#### betaMcmcPost
#' Posterior estimates of the parameter beta for each MCMC run.
#' 
#' @title posterior estimates of the parameter beta for each MCMC run
#' @param x a \code{\linkS4class{MgsaMcmcResults}}.
#' @return \code{matrix}: Posterior estimates of the parameter beta for each MCMC run.
#' @rdname betaMcmcPost-methods
#' @exportMethod betaMcmcPost
setGeneric( "betaMcmcPost", function(x) standardGeneric( "betaMcmcPost" ) )

#' @rdname betaMcmcPost-methods
setMethod(
		"betaMcmcPost",
		signature=c( "MgsaMcmcResults" ),
		function( x ) x@betaMcmcPost
)


#### pMcmcPost
#' Posterior estimates of the parameter p for each MCMC run.
#' 
#' @title posterior estimates of the parameter p for each MCMC run
#' @param x a \code{\linkS4class{MgsaMcmcResults}}.
#' @return \code{matrix}: Posterior estimates of the parameter p for each MCMC run.
#' @rdname pMcmcPost-methods
#' @exportMethod pMcmcPost
setGeneric( "pMcmcPost", function(x) standardGeneric( "pMcmcPost" ) )

#' @rdname pMcmcPost-methods
setMethod(
		"pMcmcPost",
		signature=c( "MgsaMcmcResults" ),
		function( x ) x@pMcmcPost
)

#### setsMcmcPost
#' Posterior estimates of the set marginal probabilities for each MCMC run.
#' 
#' @title posterior estimates of the the set marginal probabilities  for each MCMC run
#' @param x a \code{\linkS4class{MgsaMcmcResults}}.
#' @return \code{matrix}: Posterior estimates of the set marginal probabilities for each MCMC run.
#' @rdname setsMcmcPost-methods
#' @exportMethod setsMcmcPost
setGeneric( "setsMcmcPost", function(x) standardGeneric( "setsMcmcPost" ) )

#' @rdname setsMcmcPost-methods
setMethod(
		"setsMcmcPost",
		signature=c( "MgsaMcmcResults" ),
		function( x ) x@setsMcmcPost
)


######## show
#' Show an \code{\linkS4class{MgsaResults}}.
#'
#' @title Show an MgsaResults
#' @param object an instance of class \code{\linkS4class{MgsaResults}}.
#' @return an invisible \code{NULL}
#' @exportMethod show
setMethod(
		"show",
		signature=c( "MgsaResults" ),
		function( object ){
			cat(
					"Object of class ",
					class( object ),
					"\n",
					object@populationSize,
					" unique elements in population.\n",
					object@studySetSizeInPopulation,
					" unique elements both in study set and in population.\n",
					sep = ""
			)
			
			cat("\nPosterior on set activity (decreasing order):\n")
			nrowShow <- min (10 , nrow(object@setsResults) )
			print( dottedTable(object@setsResults[ order(object@setsResults$estimate, decreasing = TRUE)[1:nrowShow], ] , ncols=6 ) )
			if(nrowShow < nrow(object@setsResults) ){
				cat("... and ", nrow(object@setsResults) - nrowShow, " other sets.\n" )
			}
		}
)


######## plot
#' 
#' Plot method for MgsaResults objects
#' @importFrom graphics plot
#' @param x a \code{\link{MgsaResults}} 
#' @param y unused
#' @param ... unused
#' @exportMethod plot 
setMethod(
		"plot",
		signature=c( "MgsaResults" ),
		function( x, y, ... ){
			require(gplots)
			
			nrowShow <- min (10 , nrow(x@setsResults) )
			sr = x@setsResults[ rev(order(x@setsResults$estimate, decreasing = TRUE)[1:nrowShow]), ]
			
			split.screen(c(2,1))  
			split.screen(c(1,2), screen=1)
			split.screen(c(1,3), screen=2)
			
			## bar plot top ones
			screen(3)
			par(mar = c(5, 8, 4, 2) + 0.1)
			barplot2(
					sr$estimate,
					names.arg=rownames(sr),
					las = 1,
					col = "white",
					plot.ci = TRUE,
					ci.l= sr$estimate - sr$std.error,
					ci.u = sr$estimate + sr$std.error,
					horiz = TRUE,
					xlab = "Posterior (+/- std error)",
					xlim = c(0,1)
			)
			## sets 
			screen(4)
			plot( x@setsResults$estimate, xlab="Set", ylab="Posterior" )
			
			## p
			screen(5)
			with( x@pPost, plot( value, estimate, xlab="p", ylab="Posterior" ) )
			
			## alpha
			screen(6)
			with( x@alphaPost, plot( value, estimate, xlab=expression(alpha), ylab="Posterior" ) )
			
			## beta
			screen(7)
			with( x@betaPost, plot( value, estimate, xlab=expression(beta), ylab="Posterior" ) )
			
			close.screen(all=T)
		}
)

