## the class of MgsaResults: Class of the returned value of mgsa()
## Useful for plotting and show purposes

######## Class definitions
#### MgsaResults

#'
#' Note that representation is subject to change. Use the appropriate
#' accessor function to access the attributes/slots.
#' 
#' Here are some internal details of the current implementation.
#'  
#' Attributes alphaPost, betaPost, and pPost are data frames, of which the first column
#' represents a realization value, and subsequent columns the posterior probabilities of
#' each MCMC run for that realization.
#' 
#' Attribute setsResults is a dataFrame in which the first columns represent the marginal
#' posteriors of each set. The last two hold, for each set, the counts.
#' 
#' The current representation is tightly coupled to an MCMC solver, which may
#' or may not replaced in the future. Therefore, everything is subject to change.
#' You've been warned!
#' 
#' 
#'  
#' @seealso 
#' \code{\link{populationSize}}
#' \code{\link{studySetSizeInPopulation}}
#' \code{\link{alphaPost}}
#' \code{\link{betaPost}}
#' \code{\link{pPost}}
#' 

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

setGeneric( "populationSize", function(x) standardGeneric( "populationSize" ) )

#' 
#' Returns the size of the population of which the analysis was run.
#' 
setMethod(
		"populationSize",
		signature( "MgsaResults" ),
		function( x ) x@populationSize
)

### studySetSizeInPopulation

setGeneric( "studySetSizeInPopulation", function(x) standardGeneric( "studySetSizeInPopulation" ) )

#'
#' Returns the size of the study set of which the analysis was run.
#' 
setMethod(
		"studySetSizeInPopulation",
		signature( "MgsaResults" ),
		function( x ) x@studySetSizeInPopulation
)


### alphaPost

setGeneric( "alphaPost", function(x) standardGeneric( "alphaPost" ) )

setMethod(
		"alphaPost",
		signature( "MgsaResults" ),
		function( x ) x@alphaPost
)

### betaPost

setGeneric( "betaPost", function(x) standardGeneric( "betaPost" ) )

setMethod(
		"betaPost",
		signature( "MgsaResults" ),
		function( x ) x@betaPost
)

### pPost

setGeneric( "pPost", function(x) standardGeneric( "pPost" ) )

setMethod(
		"pPost",
		signature( "MgsaResults" ),
		function( x ) x@pPost
)

### setsResults

setGeneric( "setsResults", function(x) standardGeneric( "setsResults" ) )

setMethod(
		"setsResults",
		signature( "MgsaResults" ),
		function( x ) x@setsResults
)


#'
#' Instances of this class are used to hold the additional information 
#' that was provided by running an MCMC algorithm.
#' 
#' @seealso
#'  \code{\link{steps}}
#'  \code{\link{restarts}}
#'  \code{\link{alphaMcmcPost}}
#'  \code{\link{betaMcmcPost}}
#'  \code{\link{pMcmcPost}}
#' 
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
setGeneric( "steps", function(x) standardGeneric( "steps" ) )

setMethod(
		"steps",
		signature( "MgsaMcmcResults" ),
		function( x ) x@steps
)

#### Number of restarts
setGeneric( "restarts", function(x) standardGeneric( "restarts" ) )

setMethod(
		"restarts",
		signature( "MgsaMcmcResults" ),
		function( x ) x@restarts
)

#### alphaMcmcPost
setGeneric( "alphaMcmcPost", function(x) standardGeneric( "alphaMcmcPost" ) )

setMethod(
		"alphaMcmcPost",
		signature( "MgsaMcmcResults" ),
		function( x ) x@alphaMcmcPost
)


#### betaMcmcPost
setGeneric( "betaMcmcPost", function(x) standardGeneric( "betaMcmcPost" ) )

setMethod(
		"betaMcmcPost",
		signature( "MgsaMcmcResults" ),
		function( x ) x@betaMcmcPost
)


#### pMcmcPost
setGeneric( "pMcmcPost", function(x) standardGeneric( "pMcmcPost" ) )

setMethod(
		"pMcmcPost",
		signature( "MgsaMcmcResults" ),
		function( x ) x@pMcmcPost
)



######## Accessors and replacement methods
## add some only if really needed


######### print
#setMethod(
#		"print",
#		signature( "MgsaResults" ),
#		function( x, ... ){
#			cat(
#					"Object of class ",
#					class( x ),
#					"\n",
#					populationSize(x),
#					" unique elements in population.\n",
#					studySetSizeInPopulation(x),
#					" unique elements both in study set and in population.\n",
#					"The MCMC was run with ",
#					steps(x),
#					" steps.\n",
#					sep = ""
#			)
#			cat("\nPosterior on set activity:\n")
#			print(setsResults(x))
#			cat("\nPosterior on p:\n")
#			print(pPost(x))
#			cat("\nPosterior on alpha:\n")
#			print(alphaPost(x))
#			cat("\nPosterior on beta:\n")
#			print(betaPost(x))
#		}
#)


# TODO: Show mean
######## show
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

