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
	mapping<-readGAF(goa.filename)
#	load("mapping.RObj")

	# some flybase genes
#	observations<-c("vacu","vag","val","vanin-like","vap","vari","vas","vav","veg","veil","veli")

#	mgsa(getItemsIndices(mapping,observations),mapping@sets,population=getItemsIndices(mapping,observations))
	

#	getSubMapping(mapping,getItemsIndices(mapping,observations))
	
#	res<-mgsa(observations,mapping,restarts=2)
	
	
#	load("sets.RObj")
#	o<-getItemsIndices(mapping,observations)
#	sets<-mapping@sets
##	save(o,sets,file="sets.RObj")

#	print(mgsa(observations,mapping,restarts=2))

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
#	active.sets<-sample(number.of.sets,size=2)
	active.sets<-c("GO:0080090","GO:0070887")

 	hidden<-rep(0,n)
 	hidden[unlist(sets[active.sets])]<-1
	o<-hidden
 	false.positives<-runif(sum(!hidden)) < alpha
 	false.negatives<-runif(sum(hidden)) < beta
 	o[hidden][false.negatives] <- F
	o[!hidden][false.positives] <- T

	r<-mgsa(names(mapping@itemName2ItemIndex)[which(o==1)],mapping,steps=1000000)

 t<-system.time(r<-mgsa(which(o==1),mapping,steps=1000000))
 print(t)
 r

}

#mapping<-new("MgsaSets",sets=list(a=c("g1","g2"), b="g2"));
#print(mapping)

#goa.filename<-"/home/sba/.ontologizer/workspace/.cache/c5018986_0"
#mapping<-mgsa.make.go.mapping.from.goa(goa.filename)

#topgo.demo();
#mgsa.demo();
mgsa.go.demo()


sets<-list(a=c(1,2,3,4,5),d=8,e=c(2,3,4,5,6),f=c(6,7))
subset.contains<-c(1:5,8)

#
# the following juggling creates the subset mapping
#

# We assume that each set has name
if (is.null(names(sets))) names(sets)<-1:length(sets)

# First, construct a item->set mapping
# we assume that each item has at least one set
set.names<-rep(names(sets),lapply(sets,length))
set.items<-unlist(sets,use.names=F)
items<-split(set.names,set.items)

# Take the subset
items.subset<-items[subset.contains]

# Create a new set->item mapping based on the item subset
subitem.names<-rep(names(items.subset),lapply(items.subset,length))
subitem.names.f<-factor(subitem.names)
levels(subitem.names.f)<-1:length(levels(subitem.names.f))		# relabel the genes
subitem.sets<-unlist(items.subset,use.names=F)
subsets<-split(as.vector(subitem.names.f),subitem.sets)

