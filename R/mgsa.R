#'
#' Trampoline to jump into the fast C implementation.
#' 
#' @keywords internal
#' 
mgsa.trampoline <- function(o, sets, n, alpha=seq(0.01,0.3, length.out=10), beta=seq(0.1,0.8, length.out=10), p=seq(1 ,min(20,floor(length(sets)/3)), length.out=10)/length(sets), steps=1e6, restarts=1, threads=0, as=integer(0) ){
	res <- .Call("mgsa_mcmc", sets, n, o, alpha, beta, p, discrete=rep(TRUE,3), alpha.breaks=alpha, beta.breaks=beta, p.breaks=p, steps, restarts, threads, as, PACKAGE="mgsa")
	return (res)
}


#'
#' Calculates a point estimate for each row (the mean) and std.error of the mean.
#'  
#' @param x specifies a matrix of values with as a many columns as MCMC runs
#' 
#' @keywords internal
#'
mcmcSummary <- function(x){
	data.frame( estimate = rowMeans(x), std.error = apply(x,1,sd)/sqrt(ncol(x)) )	
}


#'
#' Wraps the mgsa.trampoline function. It is dumb with respect to input parameter
#' but returns a processed result.
#' 
#' @param o a vector, which defines the items that are observed. Items are expressed using
#'        integer values between 1 and \code{n}. 
#' @param sets list of sets. It is expected that each set contains integers in range between
#'        1 and \code{n} that correspond to items to which the set is associated.
#' @param n defines the total number of items. 
#' @param alpha
#' @param beta
#' @param p 
#' @param steps defines the number of MCMC steps to be performed
#' @param restarts defines the number of restarts.
#' @param threads defines number of threads to be used. Defaults to 0 which means that it
#'        corresponds to the number of available cores. 
#' @param specifies the debug level. Mainly for internal use.
#' 
#' @return an object of class MgsaMcmcResults.
#' 
#' @keywords internal 
#' 
mgsa.wrapper <- function(o, sets, n, alpha=seq(0.01,0.3, length.out=10), beta=seq(0.1,0.8, length.out=10), p=seq(1 ,min(20,floor(length(sets)/3)), length.out=10)/length(sets), steps=1e6, restarts=1, threads=0, as=integer(0), debug=0)
{
	## call to core function on non-empty sets only
	isempty  <- sapply(sets,length) == 0
	raw <- mgsa.trampoline(o, sets[!isempty], n, alpha=alpha, beta=beta, p=p, steps=steps, restarts=restarts, threads=threads, as)
	
	# just return the score (this function is quite overloaded now..., perhaps it would be better to make a separate function)
	if (length(as) > 0)
	{
		return(raw)
	}
	
	if (debug > 0)
	{
		print("Raw results:")
		str(raw);	
	}

	## wrap raw results into a MgsaResults object
	res <- new("MgsaMcmcResults")
	
	res@restarts <- restarts
	res@steps <- steps
	res@populationSize <- n
	res@studySetSizeInPopulation <- length(unique(unlist(o)))
	
	res@alphaMcmcPost <- matrix(raw$alpha$counts/steps, ncol=restarts)
	res@betaMcmcPost <- matrix(raw$beta$counts/steps, ncol=restarts)
	res@pMcmcPost <- matrix(raw$p$counts/steps, ncol=restarts)
	
	## posterior for empty sets is NA
	res@setsMcmcPost <- matrix(as.numeric(NA), nrow=length(sets), ncol=restarts)
	res@setsMcmcPost[!isempty,] <- matrix(raw$marg, ncol=restarts)
	
	res@alphaPost <- data.frame( value = raw$alpha$breaks, mcmcSummary(res@alphaMcmcPost) )
	res@betaPost <- data.frame( value = raw$beta$breaks, mcmcSummary(res@betaMcmcPost) )
	res@pPost <- data.frame( value = raw$p$breaks, mcmcSummary(res@pMcmcPost) )
	
	res@setsResults <- data.frame(
			inPopulation = sapply(sets, length),
			inStudySet = sapply(sets, function(x) length(intersect(o,x))),
			mcmcSummary(res@setsMcmcPost)
	)
	
	return(res)
}

#'
#' The main function that treats the case of character and integer inputs.
#' 
#' @param o a vector that defines the items that are observed. Items can be anything that
#'        occurs in sets.
#' @param sets a list of sets. Each set is a vector that contains all the items. The vector
#'        can be of any data type, for instance, integers or characters.
#' @param population defines the set of item that should be considered for this calculation.
#' @param alpha
#' @param beta
#' @param p 
#' @param steps defines the number of MCMC steps to be performed.
#' @param restarts defines the number of restarts.
#' @param threads defines number of threads to be used. Defaults to 0 which means that it
#'        corresponds to the number of available cores.
#' @param as define the active sets. No MCMC is run in this case. Just the log likelihood is returned.  
#'
#' @keywords internal
#' 
mgsa.main <- function(o, sets, population=NULL, alpha=seq(0.01,0.3, length.out=10), beta=seq(0.1,0.8, length.out=10), p=seq(1 ,min(20,floor(length(sets)/3)), length.out=10)/length(sets), steps=1e6, restarts=1, threads=0, as=integer(0), debug=0){
	
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

	if (debug)
	{
		cat(paste("number of sets:",length(sets),"\n"))
	}

	return(mgsa.wrapper(o,sets,length(population),alpha,beta,p,steps,restarts,threads,as,debug))
}


## S4 implementation: generic declaration
#' 
#' Performs an Mgsa analysis.
#' 
#' Currently, the Mgsa problem is solved using an MCMC sampling algorithm.  
#' 
#' @param o The observations. \code{o} can be a \code{numeric}, \code{integer}, \code{character} or \code{logical}.
#'  If \code{o} is \code{numeric}, \code{integer} or \code{character}, the entries are the items that are observed positive.
#' If \code{o} is \code{logical}, its codes the observation vector (\code{TRUE} for the items observed positives \code{FALSE} otherwise) .   
#' @param sets an instance of class MgsaSets. Alternatively, a list of sets can be specified.
#'        Each set is a vector that contains associated items. The vector can be of any data type,
#'        for instance, integers or characters.
#' @param population
#' @param alpha
#' @param beta
#' @param p
#' @param steps defines the number of the MCMC sampler. A recommended value is 1e6 or greater. 
#' @param restarts defines the number of MCMC restarts. Must be greater or equal to 1. Default to 5.
#' @param threads defines the number of threads that should be used. A value of 0 means to use all available cores. Default to 0.
#' 
#' @references GOing Bayesian: model-based gene set analysis of genome-scale data.
#' @rdname mgsa
#'  
setGeneric(
		name="mgsa",
		def=function( o, sets, population=NULL, alpha=seq(0.01,0.3, length.out=10), beta=seq(0.1,0.8, length.out=10), p=seq(1 ,min(20,floor(length(sets)/3)), length.out=10)/length(sets), steps=1e6, restarts=5, threads=0,  as=integer(0), debug=0){
			standardGeneric("mgsa")
		}
)

#' 
#' Performs an mgsa analysis
#' 
#' @param o
#' 
setMethod(
		f="mgsa",
		signature = c(o="integer", sets="list"),
		def=function( o, sets, population, alpha, beta, p, steps, restarts, threads, as, debug)
		{
			mgsa.main(o, sets, population, alpha, beta, p, steps, restarts, threads, as, debug)
		}
)

#' o numeric and sets list
setMethod(
		f="mgsa",
		signature = c(o="numeric", sets="list"),
		def=function( o, sets, population, alpha, beta, p, steps, restarts, threads, as, debug)
		{
			mgsa.main(o, sets, population, alpha, beta, p, steps, restarts, threads, as, debug)
		}
)

#' o character and sets list
setMethod(
		f="mgsa",
		signature = c(o="character", sets="list"),
		def=function( o, sets, population, alpha, beta, p, steps, restarts, threads, as, debug)
		{
			mgsa.main(o, sets, population, alpha, beta, p, steps, restarts, threads, as, debug)
		}
)

#' @rdname mgsa
#' o logical => coerce to integer with which() and call mgsa()
setMethod(
		f="mgsa",
		signature = c(o="logical", sets="list"),
		def=function( o, sets, population, alpha, beta, p, steps, restarts, threads, as, debug) {
			if (is.null(population)) population <- 1:length(o)
			mgsa( which(o), sets, population, alpha, beta, p, steps, restarts, threads, as, debug )
		}
)


setMethod(
		f="mgsa",
		signature = c(o="character", sets="MgsaSets"),
		def=function( o, sets, population, alpha, beta, p, steps, restarts, threads, as, debug) {
			if (is.null(population))
			{
				# If no population has been specified, we do not need
				# to consolidate the set and obervation ids
				return(mgsa.wrapper(getItemsIndices(sets,o), sets@sets, sets@numberOfItems, alpha, beta, p, steps, restarts, threads, as, debug))
			}
			else
			{
				population<-getItemsIndices(sets,population)
				return(mgsa ( getItemsIndices(sets,o), sets@sets, population, alpha, beta, p, steps, restarts, threads, as, debug))
			}
		}
)

#
# topGO support is optional. It is provided for convenience
#
# Hmm...any ideas for doing this in a more elegant fashion?
# (sets is not required at all here)
#if (F) # "topGO" %in% installed.packages()[,1])
#{
#	library(topGO)
#	
#	setMethod(
#			f="mgsa",
#			signature = c(o="topGOdata",sets="missing"),
#			def=function( o, sets, population, alpha, beta, p, steps, restarts, threads) {
#				data <- o
#				mgsa.main(sigGenes(data), genesInTerm(data), population, alpha, beta, p, steps, restarts, threads)
#			})
#}

