#' @include MgsaResults-class.R
#' @include MgsaSets-class.R
NULL

#' Trampoline to jump into the fast C implementation.
#' 
#' @noRd
#' @useDynLib mgsa

mgsa.trampoline <- function(o, sets, n, alpha, beta, p, steps, burnin, thin, flip.freq, restarts, threads, as ){
	res <- .Call("mgsa_mcmc", sets, n, o, alpha, beta, p, discrete=rep(TRUE,3),
				alpha.breaks=alpha, beta.breaks=beta, p.breaks=p,
				steps, burnin, thin, flip.freq,
				restarts, threads, as,
				PACKAGE="mgsa")
	return (res)
}


#'
#' Calculates a point estimate for each row (the mean) and std.error of the mean.
#'  
#' @param x specifies a matrix of values with as a many columns as MCMC runs
#' 
#' @noRd
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
#' @param p the grid defining the probability of a set to be active.
#' @param steps number of steps in each Monte-Carlo Markov chain
#' @param burnin number of burn-in MCMC steps
#' @param thin sample collecting period
#' @param restarts defines the number MCMC of restarts.
#' @param threads defines number of threads to be used. Defaults to 0 which means that it
#'        corresponds to the number of available cores.
#' @param as if not empty, the integer vector of active sets indices. No MCMC is run in this case.
#'        Just the log likelihood for a given sets assignment is returned.
#' @param debug specifies the debug level. Mainly for internal use.
#' 
#' @return an object of class \code{\linkS4class{MgsaMcmcResults}}.
#' 
#' @keywords internal 
#' @noRd
 
mgsa.wrapper <- function(o, sets, n,
					alpha=seq(0.01,0.3, length.out=10), beta=seq(0.1,0.8, length.out=10), p=seq(1 ,min(20,floor(length(sets)/3)), length.out=10)/length(sets),
					steps=1e6, burnin=0.5*steps, thin=100,
					flip.freq=0.8,
					restarts=1, threads=0, as=integer(0), debug=0)
{
	# Check parameter validity.
	if (max(alpha) > 1) stop(sprintf("Specified value %g for alpha is out of domain [0,1]",max(alpha)))
	if (min(alpha) < 0) stop(sprintf("Specified value %g for alpha is out of domain [0,1]",min(alpha)))
	if (max(beta) > 1) stop(sprintf("Specified value %g for beta is out of domain [0,1]",max(beta)))
	if (min(beta) < 0) stop(sprintf("Specified value %g for beta is out of domain [0,1]",min(beta))) 
	if (max(p) > 1) stop(sprintf("Specified value %g for p is out of domain [0,1]",max(p)))
	if (min(p) < 0) stop(sprintf("Specified value %g for p is out of domain [0,1]",min(p)))
	if (burnin > steps) stop(sprintf("Specified number of burn-in steps %d is greater than the total steps number %d",
											as.integer(burnin),as.integer(steps)))
	if (flip.freq <= 0 || flip.freq > 1) stop(sprintf("The specified frequency %g of sets state flipping Gibbs step must by in (0,1]"), flip.freq)

	## call to core function on non-empty sets only
	isempty  <- sapply(sets,length) == 0
	raw <- mgsa.trampoline(o, sets[!isempty], n, alpha=alpha, beta=beta, p=p,
					steps=steps, burnin=burnin, thin=thin, flip.freq=flip.freq,
					restarts=restarts, threads=threads, as=as)
	
	if (is.null(raw)) stop("No results could be obtained")

	# just return the score (this function is quite overloaded now..., perhaps it would be better to make a separate function)
	if (length(as) > 0)
	{
		return(raw)
	}
	
	if (debug > 0)
	{
		print("Raw results:")
		print(str(raw))	
	}

	## wrap raw results into a MgsaResults object
	res <- new("MgsaMcmcResults")
	
	res@restarts <- restarts
	res@nsamples <- raw$nsamples
	res@steps <- steps
	res@populationSize <- n
	res@studySetSizeInPopulation <- length(unique(unlist(o)))
	
	res@alphaMcmcPost <- matrix(raw$alpha$counts/raw$nsamples, ncol=restarts)
	res@betaMcmcPost <- matrix(raw$beta$counts/raw$nsamples, ncol=restarts)
	res@pMcmcPost <- matrix(raw$p$counts/raw$nsamples, ncol=restarts)
	
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
#' @param debug is used for internal debugging purposes only.
#'
#' @keywords internal
#' @noRd

mgsa.main.debug <- function(o, sets, population=NULL, debug=0, ...){
	
	if( any( sapply(sets, class)!=class(o) ) ) stop("All entries in 'sets' must have the same class as 'o'.")
	
	## population
	if(is.null(population)){
		population <- sort(unique(unlist(sets)))
		
	}else{
		if(class(population)!= class(o)){
			stop(paste("'population' must be NULL or have the same class as 'o'. The given value is",population))
		}
		population <- sort(unique(population))
	}
	
	
	## encoding (index mapping) of the elements
	encode <- function(x){ match( intersect(x, population), population) }
	sets <- lapply(sets, encode)
	no <- encode(o)

	if (length(o) != 0 && length(no) == 0) warning("None of the observations mapped to items in the sets")

	if (debug)
	{
		cat(paste("number of sets:",length(sets),"\n"))
	}

	return(mgsa.wrapper(no,sets,length(population),debug=debug,...))
}

#'
#' The main function that treats the case of character and integer inputs. This delegates
#' to mgsa.main.debug.
#' 
#' @param o a vector that defines the items that are observed. Items can be anything that
#'        occurs in sets.
#' @param sets a list of sets. Each set is a vector that contains all the items. The vector
#'        can be of any data type, for instance, integers or characters.
#' @param population defines the set of item that should be considered for this calculation.
#' @param p the grid defining the probability of a set to be active.
#'
#' @keywords internal
#' @noRd

mgsa.main <- function(o, sets, population=NULL, p=seq(1 ,min(20,floor(length(sets)/3)), length.out=10)/length(sets), ...){
	return(mgsa.main.debug(o,sets,population=population,p=p, ...))
}


## S4 implementation: generic declaration
#'  
#' Estimate marginal posterior of the MGSA problem with an MCMC sampling algorithm.  
#' 
#' The function can handle items (such as genes) encoded as \code{character} or \code{integer}.
#' For convenience \code{numeric} items can also be provided but these values should essentially be integers.
#' The type of items in the observations \code{o}, the \code{sets} and in the optional \code{population} should be consistent. 
#' In the case of \code{character} items, \code{o} and \code{population} should be of type \code{character} and \code{sets} can either be an \code{\linkS4class{MgsaSets}} or a \code{list} of \code{character} vectors.
#' In the case of \code{integer} items, \code{o} should be of type \code{integer}, \code{numeric} (but essentially with integer values), 
#' or \code{logical} and entries in \code{sets} as well as the \code{population} should be \code{integer}.
#' When \code{o} is \code{logical}, it is first coerced to integer with a call on \code{\link{which}}.
#' Observations outside the \code{population} are not taken into account. If \code{population} is \code{NULL}, it is defined as the union of all sets.
#' 
#' The default grid value for p is such that between 1 and 20 sets are active in expectation.
#' The lower limit is constrained to be lower than 0\.1 and the upper limit lower than 0\.3 independently of the total number of sets to make sure that complex solutions are penalized.
#' Marginal posteriors of activity of each set are estimated using an MCMC sampler as described in Bauer et al., 2010.
#' Because convergence of an MCM sampler is difficult to assess, it is recommended to run it several times (using \code{restarts}).
#' If variations between runs are too large (see \code{\linkS4class{MgsaResults}}), the number of steps (\code{steps}) of each MCMC run should be increased.
#' 
#' @title Performs an MGSA analysis
#' 
#' @param o The observations. It can be a \code{numeric}, \code{integer}, \code{character} or \code{logical}. See details. 
#' @param sets The sets. It can be an \code{\linkS4class{MgsaSets}} or a \code{list}. In this case, each list entry is a vector of type \code{numeric}, \code{integer}, \code{character}. See details.
#' @param population The total population. Optional. A \code{numeric}, \code{integer} or \code{character} vector.
#' Default to \code{NULL}. See details.
#' @param p Grid of values for the parameter p. Values represent probabilities of term activity and therefore must be in [0,1].
#' @param ... Optional arguments that are passed to the methods. Supported parameters are
#' \describe{
#'   \item{\code{alpha}}{Grid of values for the parameter alpha. Values represent probabilities of false-positive events and hence must be in [0,1]. \code{numeric}.}
#'   \item{\code{beta}}{Grid of values for the parameter beta. Values represent probabilities of false-negative events and hence must be in [0,1]. \code{numeric}.}
#'   \item{\code{steps}}{The number of steps of each run of the MCMC sampler. \code{integer} of length 1. A recommended value is 1e6 or greater.}
#'   \item{\code{burnin}}{The number of burn-in MCMC steps, until sample collecting begins. \code{integer} of length 1. A recommended value is half of total MCMC steps.}
#'   \item{\code{thin}}{The sample collecting period. An \code{integer} of length 1. A recommended value is 100 to reduce autocorrelation of subsequently collected samples.}
#'   \item{\code{flip.freq}}{The frequency of MCMC Gibbs step that randomly flips the state of a random set from active to inactive or vice versa. \code{numeric} from (0,1].}
#'   \item{\code{restarts}}{The number of different runs of the MCMC sampler. \code{integer} of length 1. Must be greater or equal to 1. A recommended value is 5 or greater.}
#'   \item{\code{threads}}{The number of threads that should be used for concurrent restarts. A value of 0 means to use all available cores. Default to 0.}
#' }
#' 
#' @references Bauer S., Gagneur J. and Robinson P. GOing Bayesian: model-based gene set analysis of genome-scale data. Nucleic Acids Research (2010) \url{http://nar.oxfordjournals.org/content/38/11/3523.full}
#' @return An \code{\linkS4class{MgsaMcmcResults}} object.
#' 
#' @seealso \code{\linkS4class{MgsaResults}}, \code{\linkS4class{MgsaMcmcResults}}
#' @examples
#' ## observing items A and B, with sets {A,B,C} and {B,C,D}
#' mgsa(c("A", "B"), list(set1 = LETTERS[1:3], set2 = LETTERS[2:4]))
#' 
#' ## same case with integer representation of the items and logical observation
#' mgsa(c(TRUE,TRUE,FALSE,FALSE), list(set1 = 1:3, set2 = 2:4))
#' 
#' ## a small example with gene ontology sets and plot 
#' data(example)
#' fit = mgsa(example_o, example_go)
#' ## Not run: 
#' plot(fit)
#' ## End(Not run)
#' @exportMethod mgsa
#' @rdname mgsa-methods
setGeneric(
		name="mgsa",
		def=function( o, sets, population=NULL, p=seq(min(0.1, 1/length(sets)), min(0.3, 20/length(sets)), length.out=10), ...){
			standardGeneric("mgsa")
		}
)

# o integer and sets list
#' @rdname mgsa-methods
setMethod(
		f="mgsa",
		signature = c(o="integer", sets="list"),
		def = mgsa.main
)

# o numeric and sets list
#' @rdname mgsa-methods
setMethod(
		f="mgsa",
		signature = c(o="numeric", sets="list"),
		def = mgsa.main
)

# o character and sets list
#' @rdname mgsa-methods
setMethod(
		f="mgsa",
		signature = c(o="character", sets="list"),
		def = mgsa.main
)


# o logical => coerce to integer with which() and call mgsa()
#' @rdname mgsa-methods
setMethod(
		f="mgsa",
		signature = c(o="logical", sets="list"),
		def=function( o, sets, ...) {
			mgsa( which(o), sets, ...)
		}
)

# o character and sets MgsaSets
#' @rdname mgsa-methods
setMethod(
		f="mgsa",
		signature = c(o="character", sets="MgsaSets"),
		def=function( o, sets, population=NULL, 
			      p=seq(min(0.1, 1/length(sets)), min(0.3, 20/length(sets)), length.out=10), ...) {
			
			items<-itemIndices(sets,o)
			items.nas<-sum(is.na(items))
			if (items.nas > 0)
			{
				# filter out any NAs produced during the mapping.
				# Inform the user about this procedure.
				warning(sprintf("Specified observation contains %d unmapple items. Excluded them from calculation.", items.nas))
				items<-na.omit(items)
			}

			if (is.null(population))
			{
				# If no population has been specified, we do not need
				# to consolidate the set and obervation ids
				rv = mgsa.wrapper(items, sets@sets, sets@numberOfItems, ...)
			}
			else
			{
				population<-itemIndices(sets,population)
				rv = mgsa.main( items, sets@sets, population, ...)
			}
			rv@setsResults = cbind(setsResults(rv),  setAnnotations(sets)[ rownames(setsResults(rv) ), , drop=FALSE] )
			return( rv )
		}
)
