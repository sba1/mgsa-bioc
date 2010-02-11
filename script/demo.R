#
# test with, e.g.,
#  R CMD INSTALL ../workspace/mgsa ; R --vanilla <../workspace/mgsa/script/demo.R

library(mgsa)

# TODO: Find out why functions are no longer defined
source("../workspace/mgsa/R/mgsa_make_mapping.R")

# Performs tests of random data
#
# e.g. mgsa.go.demo("/home/sba/.ontologizer/workspace/.cache/c5018986_0")
#
mgsa.go.demo<-function(goa.filename, gene.id.col = 3, go.id.col = 5, evidence.col =  7)
{
	goa.filename<-"/home/sba/.ontologizer/workspace/.cache/c5018986_0"
	mapping<-mgsa.make.go.mapping.from.goa(goa.filename)
	
	# some flybase genes
	observations<-c("vacu","vag","val","vanin-like","vap","vari","vas","vav","veg","veil","veli")
	mgsa(observations,mapping)
	
	#
	# random
	#
	
	# initialization stuff
	number.of.sets<-length(mapping@sets)
	n<-length(mapping@item.idx.map)
	sets<-mapping@sets

	alpha<-0.05
	beta<-0.05
	
	# choose two terms
	active.sets<-sample(number.of.sets,size=2)

 	hidden<-rep(0,n)
 	hidden[unlist(sets[active.sets])]<-1
	o<-hidden
 	false.positives<-runif(sum(!hidden)) < alpha
 	false.negatives<-runif(sum(hidden)) < beta
 	o[hidden][false.negatives] <- F
	o[!hidden][false.positives] <- T

 t<-system.time(r<-mgsa(which(o==1),sets,steps=1000000))
 print(t)
 r

}

mgsa.go.demo();
