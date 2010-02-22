#
# test with, e.g.,
#  R CMD INSTALL ../workspace/mgsa ; R --vanilla <../workspace/mgsa/script/demo.R

library(mgsa)


# Demonstration using topGO annotation data
topgo.demo<-function()
{
	library(topGO)
	library(ALL)
	data(ALL)
	data(geneList)
	affyLib <- paste(annotation(ALL), "db", sep = ".")
	library(package = affyLib, character.only = TRUE)
	sum(topDiffGenes(geneList))
	sampleGOdata <- new("topGOdata", description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.db,affyLib = affyLib)
	data<-sampleGOdata

	o<-sigGenes(data)
	sets<-genesInTerm(data)
	
	print(mgsa(data),restarts=1)	
}

mgsa.demo<-function()
{
	res<-mgsa(c(1,2),list(a=c(1,2),b=c(3)),steps=1e6,restarts=1)
	print(str(res))
	print(res@setsResults)
	print(alphaPost(res))
	print(betaPost(res))
	print(pPost(res))
	print(setsResults(res))
	show(res)
	
}

# Performs tests of random data
#
# e.g. mgsa.go.demo("/home/sba/.ontologizer/workspace/.cache/c5018986_0")
#
mgsa.go.demo<-function(goa.filename, gene.id.col = 3, go.id.col = 5, evidence.col =  7)
{
	# Basic
	
	### Profiling
#	Rprof("mgsa_go_demo.Rprof"); 

	goa.filename<-"/home/sba/.ontologizer/workspace/.cache/c5018986_0"
	mapping<-mgsa.make.go.mapping.from.goa(goa.filename)
#	load("mapping.RObj")

	# some flybase genes
	observations<-c("vacu","vag","val","vanin-like","vap","vari","vas","vav","veg","veil","veli")
	res<-mgsa(observations,mapping,restarts=2)
	
	
#	load("sets.RObj")
#	o<-getItemsIndices(mapping,observations)
#	sets<-mapping@sets
##	save(o,sets,file="sets.RObj")

	print(mgsa(observations,mapping,restarts=2))

	### Profiling
#	print(str(res))
#	Rprof()
#	print(summaryRprof("mgsa_go_demo.Rprof")) 
	
	#
	# random
	#
	
	# initialization stuff
	number.of.sets<-length(mapping@sets)
	n<-mapping@numberOfItems
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

#mapping<-new("MgsaMapping",sets=list(a=c("g1","g2"), b="g2"));
#print(mapping)

#goa.filename<-"/home/sba/.ontologizer/workspace/.cache/c5018986_0"
#mapping<-mgsa.make.go.mapping.from.goa(goa.filename)

#topgo.demo();
#mgsa.demo();
mgsa.go.demo()
