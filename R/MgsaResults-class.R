## the class of MgsaResults: Class of the returned value of mgsa()
## Useful for plotting and show purposes

######## Class definitions
#### MgsaResults

setClass(
        "MgsaResults",
        representation = representation(
							numberOfSteps = "numeric", populationSize = "numeric", studySetSizeInPopulation="numeric",
							alphaPost = "data.frame", betaPost ="data.frame", pPost ="data.frame",
							setsResults = "data.frame"
			)
)

######## Accessors and replacement methods

#### numberOfSteps
setGeneric( "numberOfSteps", function(x) standardGeneric( "numberOfSteps" ) )

setMethod(
        "numberOfSteps",
        signature( "MgsaResults" ),
        function( x ) x@numberOfSteps
)

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
		function( x ) x@alphaPost
)

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
		function( x ) x@betaPost
)

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
		function( x ) x@pPost
)

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
		function( x ) x@setsResults
)

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

