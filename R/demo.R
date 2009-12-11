#

# Performs tests of random data
#
# e.g. mgsa.go.demo("/home/sba/.ontologizer/workspace/.cache/c5018986_0")
#
mgsa.go.demo<-function(goa.filename, gene.id.col = 3, go.id.col = 5, evidence.col =  7)
{
	goa.filename<-"/home/sba/.ontologizer/workspace/.cache/c5018986_0"
	mapping<-mgsa.make.go.mapping.from.goa(goa.filename)
	
	# flybase genes
	observations<-c("vacu","vag","val","vanin-like","vap","vari","vas","vav","veg","veil","veli")

	mgsa(sets=mapping$sets,n=mapping$number.of.items,o=mapping$get.index(observations),steps=1000000)
	
	#
	# random
	#
	
	# initialization stuff
	number.of.sets<-length(mapping$sets)
	n<-mapping$number.of.items
	sets<-mapping$sets

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

 t<-system.time(r<-mgsa(sets,n,which(o==1),steps=1000000))
 print(t)
 r

}
