#
# Perform MGSA
#
mgsa<-function(sets,n,o,alpha=NA,beta=NA,p=NA,steps=100)
{
	r<-.Call("mgsa_mcmc",sets,n,o,alpha,beta,p,steps)
	return (r)
}

#
# Performs a basic test
#
mgsa.test<-function()
{
 n<-1000
 number.of.sets<-1000
 sets<-lapply(sample(100,size=number.of.sets,replace=T),function(x){return(sample(n,size=x,replace=F))})
 active.sets<-sample(number.of.sets,size=2)
 
 print(active.sets)
 
 alpha<-0.05
 beta<-0.05
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
