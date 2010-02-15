#'
#' Trampoline to jump into the fast C implementation.
#'
mgsa.trampoline <- function(o, sets, n, alpha=NA, beta=NA, p=NA, steps=1e6, restarts=1, threads=0 ){
			res <- .Call("mgsa_mcmc", sets, n, o, alpha, beta, p, steps, restarts, threads)
			return (res)
		}


#'
#' Wraps the mgsa.trampoline function. It is dumb with respect to input parameter
#' but returns a processed result.
#' 
#' @param o which items are observed. Items are represented by integers  
#' @param sets list of sets. It is expected that each set contains integers
#'        that correspond to item to which the set is associated
#' @param n total number of items.
#' @param alpha
#' @param beta
#' @param p 
#' @param steps number of MCMC steps to be performed
#' @param restarts defines the number of restarts.
#' @param threads defines number of threads to be used. Defaults to 0 which means that it
#'        corresponds to the number of available cores. 
#' 
#' @return an object of class MgsaResults.
#' 

mgsa.wrapper <- function(o, sets, n, alpha=NA, beta=NA, p=NA, steps=1e6, restarts=1, threads=0){
	## call to core function
	raw <- mgsa.trampoline(o, sets, n, alpha=alpha, beta=beta, p=p, steps=steps, restarts=restarts, threads=threads)
	
	## wrap raw results into a MgsaResults object
	res <- new("MgsaResults")
	
	res@numOfRestarts <- as.integer(restarts)
	numberOfSteps(res) <- steps
	populationSize(res) <- n
	studySetSizeInPopulation(res) <- length(unique(unlist(o)))
	
	alphaPost(res) <- data.frame( value = raw$alpha$breaks, posterior=raw$alpha$counts/steps )
	betaPost(res) <- data.frame( value = raw$beta$breaks, posterior=raw$beta$counts/steps )
	pPost(res) <- data.frame( value = raw$p$breaks, posterior = raw$p$counts/steps )
	
	setsResults(res) <- data.frame(
			posterior=raw$marg,
			inPopulation = sapply(sets, length),
			inStudySet = sapply(sets, function(x) length(intersect(o,x)))
	)
	
	return(res)
}

#'
#' The main function that treats the case of character and integer inputs.
#' 
#' @param o 
#' @param sets
#' @param population
#' @param alpha
#' @param beta
#' @param p 
#' @param steps number of MCMC steps to be performed
#' @param restarts defines the number of restarts.
#' @param threads defines number of threads to be used. Defaults to 0 which means that it
#'        corresponds to the number of available cores. 
#' 
mgsa.main <- function(o, sets, population=NULL, alpha=NA, beta=NA, p=NA, steps=1e6, restarts=1, threads=0){
	
	if( any( sapply(sets, class)!=class(o) ) ) stop("All entries in 'sets' must have the same class as 'o'.")
	
	## population
	if(is.null(population)){
		population <- sort(unique(unlist(sets)))
		
	}else{
		if(class(population)!= class(o)){
			stop("'population' must be NULL or have the same class as 'o'.")
		}
		population <- sort(unique(population))
	}
	
	## encoding (index mapping) of the elements
	encode <- function(x){ match( intersect(x, population), population) }
	sets <- lapply(sets, encode)
	o <- encode(o)

	return(mgsa.wrapper(o,sets,length(population),alpha,beta,p,steps,restarts,threads))
}


## S4 implementation: generic declaration
setGeneric(
		"mgsa",
		function( o, sets, population=NULL, alpha=NA, beta=NA, p=NA, steps=1e6, restarts=1, threads=0 ){
			standardGeneric("mgsa")
		}
)

# o integer and sets list
setMethod(
		"mgsa",
		signature = c(o="integer", sets="list"),
		function( o, sets, population, alpha, beta, p, steps, restarts, threads) mgsa.main(o, sets, population, alpha, beta, p, steps, restarts, threads)
)

# o numeric and sets list
setMethod(
		"mgsa",
		signature = c(o="numeric", sets="list"),
		function( o, sets, population, alpha, beta, p, steps, restarts, threads)
		{
			mgsa.main(o, sets, population, alpha, beta, p, steps, restarts, threads)
		}
)

# o character and sets list
setMethod(
		"mgsa",
		signature = c(o="character", sets="list"),
		function( o, sets, population, alpha, beta, p, steps, restarts, threads) mgsa.main(o, sets, population, alpha, beta, p, steps, restarts, threads)
)

# o logical => coerce to integer with which() and call mgsa()
setMethod(
		"mgsa",
		signature = c(o="logical", sets="list"),
		function( o, sets, population, alpha, beta, p, steps, restarts, threads) {
			if (is.null(population)) population <- 1:length(o)
			mgsa( which(o), sets, population, alpha, beta, p, steps, restarts, threads )
		}
)

# o character and mapping object
# TODO: make it also work with integers
setMethod(
		"mgsa",
		signature = c(o="character", sets="MgsaMapping"),
		function( o, sets, population, alpha, beta, p, steps, restarts, threads) {
			if (is.null(population))
			{
				# If no population has been specified, we do not need
				# to consolidate the set and obervation ids
				return(mgsa.wrapper( sets@item.idx.map[o], sets@sets, length(sets@item.idx.map), alpha, beta, p, steps, restarts, threads))
			}
			else
			{
				population<-sets@item.idx.map[population]
				return(mgsa ( sets@item.idx.map[o], sets@sets, population, alpha, beta, p, steps, restarts, threads))
			}
		}
)
