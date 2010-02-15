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
#' The current represenation is tightly coupled to an MCMC solver, which may
#' or may not replaced in the furture. Therefore, everything is subject to change.
#' You've been warned!
#' 

setClass(
        "MgsaResults",
        representation = representation(
							numberOfSteps = "numeric", populationSize = "numeric", studySetSizeInPopulation="numeric",
							numOfRestarts = "integer",
# TODO: Wrap these into own classes
							alphaPost = "data.frame", betaPost ="data.frame", pPost ="data.frame",
							setsResults = "data.frame"
			)
)

######## Accessors and replacement methods

#### numberOfSteps, MCMC specific
setGeneric( "numberOfSteps", function(x) standardGeneric( "numberOfSteps" ) )

setMethod(
        "numberOfSteps",
        signature( "MgsaResults" ),
        function( x ) x@numberOfSteps
)


# Do we need this? This is a convenience function, but a user should not be able
# to set this attribute.
setGeneric( "numberOfSteps<-", function( x, value ) standardGeneric( "numberOfSteps<-" ) )

setReplaceMethod(
        "numberOfSteps", "MgsaResults",
        function( x, value ) {
            x@numberOfSteps <- value
            return(x)
        }
)

#### populationSize
setGeneric( "populationSize", function(x) standardGeneric( "populationSize" ) )

setMethod(
		"populationSize",
		signature( "MgsaResults" ),
		function( x ) x@populationSize
)

# FIXME: Do we need this? This is a convenience function, but a user should not be able
# to set this attribute.
setGeneric( "populationSize<-", function( x, value ) standardGeneric( "populationSize<-" ) )

setReplaceMethod(
		"populationSize", "MgsaResults",
		function( x, value ) {
			x@populationSize <- value
			return(x)
		}
)



#### studySetSizeInPopulation
setGeneric( "studySetSizeInPopulation", function(x) standardGeneric( "studySetSizeInPopulation" ) )

setMethod(
		"studySetSizeInPopulation",
		signature( "MgsaResults" ),
		function( x ) x@studySetSizeInPopulation
)

# FIXME: Do we need this? This is a convenience function, but a user should not be able
# to set this attribute.
setGeneric( "studySetSizeInPopulation<-", function( x, value ) standardGeneric( "studySetSizeInPopulation<-" ) )

setReplaceMethod(
		"studySetSizeInPopulation", "MgsaResults",
		function( x, value ) {
			x@studySetSizeInPopulation <- value
			return(x)
		}
)



#### alphaPost
setGeneric( "alphaPost", function(x) standardGeneric( "alphaPost" ) )

setMethod(
		"alphaPost",
		signature( "MgsaResults" ),
		function( x ) data.frame(value=x@alphaPost[,1],posterior=rowMeans(x@alphaPost[,-1]))
)

# FIXME: Do we need this? This is a convenience function, but a user should not be able
# to set this attribute.
setGeneric( "alphaPost<-", function( x, value ) standardGeneric( "alphaPost<-" ) )

setReplaceMethod(
		"alphaPost", "MgsaResults",
		function( x, value ) {
			x@alphaPost <- value
			return(x)
		}
)


#### betaPost
setGeneric( "betaPost", function(x) standardGeneric( "betaPost" ) )

setMethod(
		"betaPost",
		signature( "MgsaResults" ),
		function( x ) data.frame(value=x@betaPost[,1],posterior=rowMeans(x@betaPost[,-1]))
)

# FIXME: Do we need this? This is a convenience function, but a user should not be able
# to set this attribute.
setGeneric( "betaPost<-", function( x, value ) standardGeneric( "betaPost<-" ) )

setReplaceMethod(
		"betaPost", "MgsaResults",
		function( x, value ) {
			x@betaPost <- value
			return(x)
		}
)


#### pPost
setGeneric( "pPost", function(x) standardGeneric( "pPost" ) )

setMethod(
		"pPost",
		signature( "MgsaResults" ),
		function( x ) data.frame(value=x@pPost[,1],posterior=rowMeans(x@pPost[,-1]))
)

# FIXME: Do we need this? This is a convenience function, but a user should not be able
# to set this attribute.
setGeneric( "pPost<-", function( x, value ) standardGeneric( "pPost<-" ) )

setReplaceMethod(
		"pPost", "MgsaResults",
		function( x, value ) {
			x@pPost <- value
			return(x)
		}
)


#### setsResult
setGeneric( "setsResults", function(x) standardGeneric( "setsResults" ) )

setMethod(
		"setsResults",
		signature( "MgsaResults" ),
		function( x ) data.frame(x@setsResults[,c(x@numOfRestarts+1,x@numOfRestarts+2)],posterior=rowMeans(x@setsResults[,1:x@numOfRestarts]))
)

# FIXME: Do we need this? This is a convenience function, but a user should not be able
# to set this attribute.
setGeneric( "setsResults<-", function( x, value ) standardGeneric( "setsResults<-" ) )

setReplaceMethod(
		"setsResults", "MgsaResults",
		function( x, value ) {
			x@setsResults <- value
			return(x)
		}
)



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
#					numberOfSteps(x),
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
					populationSize(object),
					" unique elements in population.\n",
					studySetSizeInPopulation(object),
					" unique elements both in study set and in population.\n",
					"The MCMC was run with ",
					numberOfSteps(object),
					" steps.\n",
					sep = ""
			)
			cat("\nPosterior on set activity (decreasing order):\n")
			nrowShow <- min (10 , nrow(setsResults(object)) )
			print( setsResults(object)[ order(setsResults(object)$posterior, decreasing = TRUE)[1:nrowShow], ] )
			if(nrowShow < nrow(setsResults(object)) ){
				cat("... and ", nrow(setsResults(object)) - nrowShow, " other sets.\n" )
			}
		}
)

# TODO: Plot mean
######## plot
setMethod(
		"plot",
		signature( "MgsaResults" ),
		function( x, y, ... ){
			par( mfrow=c(2,2) )
			## sets 
			plot( setsResults(x)$posterior, xlab="Set", ylab="Posterior" )
			## p 
			with( pPost(x), plot( value, posterior, xlab="p", ylab="Posterior" ) )
			## alpha 
			with( alphaPost(x), plot( value, posterior, xlab=expression(alpha), ylab="Posterior" ) )
			## beta 
			with( betaPost(x), plot( value, posterior, xlab=expression(beta), ylab="Posterior" ) )
		}
)

