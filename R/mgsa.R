#
# Perform MGSA
#
mgsa.core <- function(o, sets, n, alpha=NA, beta=NA, p=NA, steps=1e6 ){
			res <- .Call("mgsa_mcmc", sets, n, o, alpha, beta, p, steps)
			return (res)
		}


## main function that treats case of character and integer
mgsa.main <- function(o, sets, population=NULL, alpha=NA, beta=NA, p=NA, steps=1e6){
	
	if( any( sapply(sets, class)!=class(o) ) ) stop("All entries in 'sets' must have the same class than 'o'.")
	
	## population
	if(is.null(population)){
		population <- sort(unique(unlist(sets)))
		
	}else{
		if(class(population)!= class(o)){
			stop("'population' must be NULL or have the same class than 'o'.")
		}
		population <- sort(unique(population))
	}
	
	encode <- function(x){ match( intersect(x, population), population) }
	sets <- lapply(sets, encode)
	o <- encode(o)
	res <- mgsa.core(o, sets, n=length(population), alpha=alpha, beta=beta, p=p, steps=steps)
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


