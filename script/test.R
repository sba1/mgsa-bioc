#
# Performs a basic test
#
###############################################################################

require(mgsa)

n  <-  1000
number.of.sets <- 100
number.of.steps <- 1e6

sets <- lapply(
		sample(100,size=number.of.sets,replace=T),
		function(x){
			return(sample(n,size=x,replace=F))
		}
)

active.sets <- sample(number.of.sets,size=2)

print(active.sets)

alpha <- 0.05
beta <- 0.05
hidden <- rep(FALSE,n)
hidden[unlist(sets[active.sets])] <- TRUE
o <- hidden
false.positives <- runif(sum(!hidden)) < alpha
false.negatives <- runif(sum(hidden)) < beta
o[hidden][false.negatives]  <-  FALSE
o[!hidden][false.positives]  <-  TRUE

## integer, list
cat("mgsa: integer, list:\n")
#t <- system.time(r <- mgsa(which(o==1), sets, 1:n, steps=1e6))
t <- system.time(r <- mgsa(which(o==1), sets, steps=number.of.steps))
print(t)
print(r)
plot(r)


## integer, list, setting alpha and beta grids
## Not yet implemented
cat("mgsa: integer, list, setting alpha and beta grids:\n")
r <- mgsa(which(o==1), sets, steps=number.of.steps, alpha = 0.01*1:10, beta = 0.01*1:10)
print(r)
plot(r)


## from now on with set names
names(sets) <- paste("set",1:length(sets), sep="_")
cat("mgsa: integer, list:\n")
#t <- system.time(r <- mgsa(which(o==1), sets, 1:n, steps=1e6))
t <- system.time(r <- mgsa(which(o==1), sets, steps=number.of.steps))
print(t)
print(r)
plot(r)


## with gene names
genes = sapply( 1:n, function(i) do.call( paste, c( as.list( sample( LETTERS, 6, replace=TRUE) ), sep="" ) ) )

sets2 = lapply(sets, function(x) genes[x] )
o2 = genes[o]

cat("mgsa: character, list:\n")
#t <- system.time(r <- mgsa(which(o==1), sets, 1:n, steps=number.of.steps))
t <- system.time(r <- mgsa(o2, sets2, steps=number.of.steps))
print(t)
print(r)
plot(r)

cat("mgsa: logical, list:\n")
#t <- system.time(r <- mgsa(which(o==1), sets, 1:n, steps=1e6))
t <- system.time(r <- mgsa(o==1, sets, steps=number.of.steps))
print(t)
print(r)
plot(r)

## several 
## on purpose small length to be able to see the error bars
r <- mgsa(o==1, sets, steps=1e3, restarts=10)
print(r)
plot(r)

## integer, list with a few empty sets
cat("mgsa: integer, list, with an empty set in second position:\n")
r <- mgsa(which(o==1), c(sets[1], list(integer(0)), sets[-1]), steps=1e3, restarts=10)
print(r)
plot(r)
print(head(r@setsResults))



