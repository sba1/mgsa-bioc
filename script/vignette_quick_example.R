# create example for the vignette

DATADIR = "/localdata/mgsa"

## paper example
library(mgsa)

go = readGAF(file.path (DATADIR,"gene_association.sgd"))
number.of.sets = 10
set.seed(0)
rnd =  sample(length(go@sets), number.of.sets)
pop = unique(unlist(go@sets[rnd]))
n = length(pop)
subgo = getSubMapping(go, pop)

subgo@sets = subgo@sets[ names(subgo@sets) %in% names(go@sets)[rnd] ]
subgo@setAnnotations = go@setAnnotations[rnd, ,drop=FALSE]
subgo@itemAnnotations = go@itemAnnotations[pop, ,drop=FALSE]
names(subgo@itemName2ItemIndex) = names(go@itemName2ItemIndex)[pop]
		
## active sets
active.sets <- c(2,6)
alpha <- 0.1
beta <- 0.20
hidden <- rep(FALSE,n)
hidden[unlist(subgo@sets[active.sets])] <- TRUE
o <- hidden
false.positives <- runif(sum(!hidden)) < alpha
false.negatives <- runif(sum(hidden)) < beta
o[hidden][false.negatives]  <-  FALSE
o[!hidden][false.positives]  <-  TRUE


example_go = subgo
example_o = names(subgo@itemName2ItemIndex)[o]

## save example dataset
save(example_go, example_o, file="inst/data/example.rda")


