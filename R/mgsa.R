#
# Perform MGSA
#
mgsa<-function(sets,n,o,alpha=NA,beta=NA,p=NA,steps=100)
{
	r<-.Call("mgsa_mcmc",sets,n,o,alpha,beta,p,steps)
	return (r)
}
