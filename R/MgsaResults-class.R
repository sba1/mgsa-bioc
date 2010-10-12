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
#' The columns of the slot \code{setsResults} contains the number of items of the set in the population, the numebr of items of the set in the study set, the estimate of its marginal posterior probability and its standard error.
#' The \code{\link{rownames}} are the names of the sets if available.
#' 
#' Accessor methods exist for each slot.
#' 
#' @title Results of an MGSA analysis
#' 
#' @seealso mgsa
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

#' @exportMethod populationSize
setGeneric( "populationSize", function(x) standardGeneric( "populationSize" ) )

#' The size of the population of which the analysis was run.
#' @param x an MgsaResults object
setMethod(
		"populationSize",
		signature( "MgsaResults" ),
		function( x ) x@populationSize
)

### studySetSizeInPopulation
#' @exportMethod studySetSizeInPopulation
setGeneric( "studySetSizeInPopulation", function(x) standardGeneric( "studySetSizeInPopulation" ) )

#' The size of the study set of which the analysis was run.
setMethod(
		"studySetSizeInPopulation",
		signature( "MgsaResults" ),
		function( x ) x@studySetSizeInPopulation
)


### alphaPost
#' @exportMethod alphaPost
setGeneric( "alphaPost", function(x) standardGeneric( "alphaPost" ) )

setMethod(
		"alphaPost",
		signature( "MgsaResults" ),
		function( x ) x@alphaPost
)

### betaPost
#' @exportMethod betaPost
setGeneric( "betaPost", function(x) standardGeneric( "betaPost" ) )

setMethod(
		"betaPost",
		signature( "MgsaResults" ),
		function( x ) x@betaPost
)

### pPost
#' @exportMethod pPost
setGeneric( "pPost", function(x) standardGeneric( "pPost" ) )

setMethod(
		"pPost",
		signature( "MgsaResults" ),
		function( x ) x@pPost
)

### setsResults
#' @exportMethod setsResults
setGeneric( "setsResults", function(x) standardGeneric( "setsResults" ) )

#' @rdname MgsaResults
setMethod(
		"setsResults",
		signature( "MgsaResults" ),
		function( x ) x@setsResults
)


#'
#' Instances of this class are used to hold the additional information 
#' that was provided by running an MCMC algorithm.
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
#' @seealso mgsa
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
#' @exportMethod steps
#' @rdname MgsaResults
setGeneric( "steps", function(x) standardGeneric( "steps" ) )

setMethod(
		"steps",
		signature( "MgsaMcmcResults" ),
		function( x ) x@steps
)

#### Number of restarts
#' @exportMethod restarts
setGeneric( "restarts", function(x) standardGeneric( "restarts" ) )

#' @rdname MgsaResults
setMethod(
		"restarts",
		signature( "MgsaMcmcResults" ),
		function( x ) x@restarts
)

#### alphaMcmcPost
#' @exportMethod alphaMcmcPost
setGeneric( "alphaMcmcPost", function(x) standardGeneric( "alphaMcmcPost" ) )

#' @rdname MgsaResults
setMethod(
		"alphaMcmcPost",
		signature( "MgsaMcmcResults" ),
		function( x ) x@alphaMcmcPost
)


#### betaMcmcPost
#' @exportMethod betaMcmcPost
setGeneric( "betaMcmcPost", function(x) standardGeneric( "betaMcmcPost" ) )

#' @rdname MgsaResults
setMethod(
		"betaMcmcPost",
		signature( "MgsaMcmcResults" ),
		function( x ) x@betaMcmcPost
)


#### pMcmcPost
#' @exportMethod pMcmcPost
setGeneric( "pMcmcPost", function(x) standardGeneric( "pMcmcPost" ) )

#' @rdname MgsaResults
setMethod(
		"pMcmcPost",
		signature( "MgsaMcmcResults" ),
		function( x ) x@pMcmcPost
)

#### setsMcmcPost
#' @exportMethod setsMcmcPost
setGeneric( "setsMcmcPost", function(x) standardGeneric( "setsMcmcPost" ) )

#' @rdname MgsaResults
setMethod(
		"setsMcmcPost",
		signature( "MgsaMcmcResults" ),
		function( x ) x@setsMcmcPost
)


######## show
#' 
#' Show method for MgsaResults objects
#' @param object a \code{\link{MgsaResults}} 
#' @exportMethod show 
#' @rdname MgsaResults
setMethod(
		"show",
		signature( "MgsaResults" ),
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
			
			print(str(object@setsResults))
			
			cat("\nPosterior on set activity (decreasing order):\n")
			nrowShow <- min (10 , nrow(object@setsResults) )
			print( object@setsResults[ order(object@setsResults$estimate, decreasing = TRUE)[1:nrowShow], ] )
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
#' @exportMethod plot 
#' @rdname MgsaResults
setMethod(
		"plot",
		signature( "MgsaResults" ),
		function( x, y, ... ){
			require(gplots)
			#par( mfrow=c(2,3) )
			
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
		}
)

