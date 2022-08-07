## Code CVM 
CVM <- function(A,nbgrille,sigma)
{
	# A : sample of data
	# nbgrille : preision for the approximation of the grid
	# variance to be used in the normal distribution function
	
	## Computes sqrt(n) * \integral( F_n(x) - Phi(x) )^2 dx where Phi is the normal distribution function mean 0
	## variance sigma, F_n is the empirical distribution function computed from A
	cv =0
	n <- length(A)
		for(i in 1:nbgrille){
		cv <- cv +  n/nbgrille * (mean(A <= qnorm(i/nbgrille,mean=0,sd=sqrt(sigma)) )- i/nbgrille )^2 
	}
	cv
}
pval <- function(dis,cv){
	mean(cv>dis)
}


#Linear fit with CVM statistic
fit1CVM <- function(A,poids,b,nbgrille){
don <- as.matrix(A)
nbr <- ncol(A)-1
nbniv <- nrow(don)	
	
# Weights
W <- diag(rep(1/abs(don[,1])^poids,nbr))	
# Optimal parameters
x1 <- don[,1]
Xp <- matrix(c(rep(1,length(don[,1])*nbr),rep(x1,nbr)),ncol=2)
matInv <- solve(t(Xp)%*%W%*%Xp)
param.optimaux <- rev(matInv%*%t(Xp)%*%W%*%c(don[,-1]))

#Residuals
temp <- param.optimaux[1]*don[,1]+  param.optimaux[2]
temp2 <- matrix(rep(temp,nbr),nrow=length(don[,1]))
temp3 <- (don[,-1] - temp2)
	
poidsob <- 1/abs(don[,1])^poids
poidsob <- poidsob/sum(poidsob)
epsilonp <- c(temp3)*(sqrt(rep(poidsob,nbr)))
	
var.param <- matInv * var(epsilonp)	
	
## Bootstrap
	n <- length(epsilonp)
	sigma <- sqrt(var(epsilonp))
	distrCv <- rep(0,b)
	
for(j in 1:b){
	
	B <- rnorm(n,0,sqrt(sigma))
	C <- rmvnorm(1,mean=c(0,0),var.param)
	drift <- 	(C[1] + C[2]*mean(don[,1])) # Drift function to used when esimated residuals are involved
	cv <- 0
	for(i in 1:nbgrille){
		cv <- cv +  n/nbgrille * (mean(B <= qnorm(i/nbgrille,mean=0,sd=sqrt(sigma)) )- i/nbgrille + drift*dnorm(qnorm(i/nbgrille,0,sqrt(sigma)),0,sd=sqrt(sigma))/sqrt(n) )^2 
	}
	distrCv[j] <- cv	
	}	
cv0 <- CVM(epsilonp,nbgrille,sigma)	
	pvalN <- pval(cv0,distrCv)
	
liste <- list(param=param.optimaux,varparam=var.param,pvalN=pvalN)
return(liste)
}	

# Quadratic fit statistic CVM
fit2CVM <- function(A,poids,b,nbgrille){
don <- as.matrix(A)
nbr <- ncol(A)-1
nbniv <- nrow(A)	


#Weights
W <- diag(rep(1/abs(don[,1])^poids,nbr))	
#Optimal paramters
x1 <- don[,1]
Xp <- matrix(c(rep(1,length(don[,1])*nbr),rep(x1,nbr),rep(x1^2,nbr)),ncol=3)
matInv <- solve(t(Xp)%*%W%*%Xp)
param.optimaux <- rev(matInv%*%t(Xp)%*%W%*%c(don[,-1]))

	
#Residuals
temp <- param.optimaux[1]*don[,1]^2 +  param.optimaux[2]*don[,1] +  param.optimaux[3]
temp2 <- matrix(rep(temp,nbr),nrow=nbniv)
temp3 <- (don[,-1] - temp2)

poidsob <- 1/abs(don[,1])^poids
poidsob <- poidsob/sum(poidsob)
epsilonp <- c(temp3)*(sqrt(rep(poidsob,nbr)))	
var.param <- matInv * var(epsilonp)
	
## Bootstrap
	n <- length(epsilonp)
	sigma <- sqrt(var(epsilonp))
	distrCv <- rep(0,b)
	
for(j in 1:b){
	
	B <- rnorm(n,0,sqrt(sigma))
	C <- rmvnorm(1,mean=c(0,0,0),var.param)
	drift <- 	(C[1] + C[2]*mean(don[,1]) + C[3]*mean(don[,1]^2)) # Drift function to used when esimated residuals are involved
	cv <- 0
	for(i in 1:nbgrille){
		cv <- cv +  n/nbgrille * (mean(B <= qnorm(i/nbgrille,mean=0,sd=sqrt(sigma)) )- i/nbgrille + drift*dnorm(qnorm(i/nbgrille,0,sqrt(sigma)),0,sd=sqrt(sigma))/sqrt(n) )^2 
	}
	distrCv[j] <- cv	
	}	
cv0 <- CVM(epsilonp,nbgrille,sigma)	
	pvalN <- pval(cv0,distrCv)
	
liste <- list(param=param.optimaux,varparam=var.param,pvalN=pvalN)
return(liste)
}	