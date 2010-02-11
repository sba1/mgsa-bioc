#
# Perform MGSA
#
mgsa.core <- function(o, sets, n, alpha=NA, beta=NA, p=NA, steps=1e6 ){
			res <- .Call("mgsa_mcmc", sets, n, o, alpha, beta, p, steps)
			return (res)
		}


## main function that treats case of character and integer
mgsa.main <- function(o, sets, population=NULL, alpha=NA, beta=NA, p=NA, steps=1e6){
	
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
	
	## call to core function
	raw <- mgsa.core(o, sets, n=length(population), alpha=alpha, beta=beta, p=p, steps=steps)
	
	## wrap raw results into a MgsaResults object
	res <- new("MgsaResults")
	
	numberOfSteps(res) <- steps
	populationSize(res) <- length(population)
	studySetSizeInPopulation(res) <- length(unique(unlist(o)))
	
	alphaPost(res) <- data.frame( value = raw$alpha$breaks, posterior = raw$alpha$counts/steps )
	betaPost(res) <- data.frame( value = raw$beta$breaks, posterior = raw$beta$counts/steps )
	pPost(res) <- data.frame( value = raw$p$breaks, posterior = raw$p$counts/steps )
	
	setsResults(res) <- data.frame(
							posterior = raw$marg,
							inPopulation = sapply(sets, length),
							inStudySet = sapply(sets, function(x) length(intersect(o,x)))
					)
	rownames( setsResults(res) ) <- names(sets)
	
	return(res)
	
}


## S4 implementation: generic declaration
setGeneric(
		"mgsa",
		function( o, sets, population=NULL, alpha=NA, beta=NA, p=NA, steps=1e6 ){
			standardGeneric("mgsa")
		}
)

# o integer and sets list
setMethod(
		"mgsa",
		signature = c(o="integer", sets="list"),
		function( o, sets, population, alpha, beta, p, steps) mgsa.main(o, sets, population, alpha, beta, p, steps)
)

# o character and sets list
setMethod(
		"mgsa",
		signature = c(o="character", sets="list"),
		function( o, sets, population, alpha, beta, p, steps) mgsa.main(o, sets, population, alpha, beta, p, steps)
)

# o logical => coerce to integer with which() and call mgsa()
setMethod(
		"mgsa",
		signature = c(o="logical", sets="list"),
		function( o, sets, population, alpha, beta, p, steps) {
			if(is.null(population)) population <- 1:length(o)
			mgsa( which(o), sets, population, alpha, beta, p, steps )
		}
)

# o character and mapping object
# TODO: make it also work with integers
setMethod(
		"mgsa",
		signature = c(o="character", sets="MgsaMapping"),
		function( o, sets, population, alpha, beta, p, steps) {
			if (!is.null(population)) population<-sets@item.idx.map[population]
			return(mgsa ( sets@item.idx.map[o], sets@sets, population, alpha, beta, p, steps))
		}
)
mgsa(observations,mapping)
