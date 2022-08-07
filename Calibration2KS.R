## Code KS

KS <- function(A,nbgrille,sigma){
	# A : sample of data
	# nbgrille : preision for the approximation of the grid
	# variance to be used in the normal distribution function
	
	## Computes sqrt(n) * sup_x | F_n(x) - Phi(x) | where Phi is the normal distribution function mean 0
	## variance sigma, F_n is the empirical distribution function computed from A
	
	n<- length(A)
	fn <- function(t){sqrt(n) * abs(mean(A <= qnorm(t,0,sd=sqrt(sigma))) - t ) }
	max(sapply((1:nbgrille)/nbgrille,fn))
}
pval <- function(dis,cv){
	mean(cv>dis)
}


# Linear fit using KS statistic
fit1KS <- function(A,poids,b,nbgrille){
don <- as.matrix(A)
nbr <- ncol(A)-1

W <- diag(rep(1/abs(don[,1])^poids,nbr))	# matrix of weights
	
# Finding the optimal parameters for the regression
x1 <- don[,1]
Xp <- matrix(c(rep(1,length(don[,1])*nbr),rep(x1,nbr)),ncol=2)
matInv <- solve(t(Xp)%*%W%*%Xp)
param.optimaux <- rev(matInv%*%t(Xp)%*%W%*%c(don[,-1]))

# Computing the residuals
temp <- param.optimaux[1]*don[,1]+  param.optimaux[2]
temp2 <- matrix(rep(temp,nbr),nrow=length(don[,1]))
temp3 <- (don[,-1] - temp2)
	
poidsob <- 1/abs(don[,1])^poids
poidsob <- poidsob/sum(poidsob)
epsilonp <- c(temp3)*(sqrt(rep(poidsob,nbr))) # residuals
	
var.param <- matInv * var(epsilonp)	# the variance of the estimated paramters
	
## Bootstrap
	n <- length(epsilonp)
	sigma <- sqrt(var(epsilonp))
	distrCv <- rep(0,b)
	
for(j in 1:b){
	
	B <- rnorm(n,0,sqrt(sigma))
	C <- rmvnorm(1,mean=c(0,0),var.param)
	drift <- 	(C[1] + C[2]*mean(don[,1])) # Drift function to used when esimated residuals are involved
	n<- length(A)
	fn <- function(t){
		sqrt(n) * abs(mean(A <= qnorm(t,0,sd=sqrt(sigma))) - 
				t +  drift*dnorm(qnorm(t,0,sqrt(sigma)),0,sd=sqrt(sigma))/sqrt(n)) }
	distrCv[j] <- max(sapply((1:nbgrille)/nbgrille,fn))	
	}	
cv0 <- KS(epsilonp,nbgrille,sigma)	
	pvalN <- pval(cv0,distrCv)
	
liste <- list(param=param.optimaux,res1=cbind(rep(poidsob,nbr),c(temp3)),epsp=epsilonp,varparam=var.param,pvalN=pvalN)
return(liste)
}	

#quadratic fit using KS statistic
fit2KS <- function(A,poids,b,nbgrille){
don <- as.matrix(A)
nbr <- ncol(A)-1

W <- diag(rep(1/abs(don[,1])^poids,nbr))	# Weights

# Optimal parameters
x1 <- don[,1]
Xp <- matrix(c(rep(1,length(don[,1])*nbr),rep(x1,nbr),rep(x1^2,nbr)),ncol=3)
matInv <- solve(t(Xp)%*%W%*%Xp)
param.optimaux <- rev(matInv%*%t(Xp)%*%W%*%c(don[,-1]))
	
#Residuals
temp <- param.optimaux[1]*don[,1]^2 +  param.optimaux[2]*don[,1] +  param.optimaux[3]
temp2 <- matrix(rep(temp,nbr),nrow=nrow(don))
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
	n<- length(A)
	fn <- function(t){
		sqrt(n) * abs(mean(A <= qnorm(t,0,sd=sqrt(sigma))) - 
				t +  drift*dnorm(qnorm(t,0,sqrt(sigma)),0,sd=sqrt(sigma))/sqrt(n)) }
	distrCv[j] <- max(sapply((1:nbgrille)/nbgrille,fn))	
	}	
cv0 <- KS(epsilonp,nbgrille,sigma)	
	pvalN <- pval(cv0,distrCv)
	pvalN <- pval(cv0,distrCv)
	
liste <- list(param=param.optimaux,varparam=var.param,pvalN=pvalN)
return(liste)
}