library(mgsa)

n  <-  1000
number.of.sets <- 100
number.of.steps <- 1e5

sets <- lapply(
		sample(100,size=number.of.sets,replace=T),
		function(x){
			return(sample(n,size=x,replace=F))
		}
)

active.sets <- sample(number.of.sets,size=2)

print(active.sets)

# generate noisy observations
alpha <- 0.05
beta <- 0.05
hidden <- rep(FALSE,n)
hidden[unlist(sets[active.sets])] <- TRUE
o <- hidden
false.positives <- runif(sum(!hidden)) < alpha
false.negatives <- runif(sum(hidden)) < beta
o[hidden][false.negatives]  <-  FALSE
o[!hidden][false.positives]  <-  TRUE


mgsa.trampoline.2 <- function(o, sets, n, alpha=NA, beta=NA, p=NA, discrete=c(T,T,T), steps=1e6, restarts=1, threads=0, as=integer(0) ){
	res <- .Call("mgsa_mcmc", sets, n, o, alpha, beta, p, discrete, steps, restarts, threads, as)
	return (res)
}

r<-mgsa.trampoline.2(which(o==1), sets, length(sets), steps=number.of.steps)
