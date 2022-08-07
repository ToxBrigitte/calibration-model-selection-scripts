	# A : sample of data
	# nbgrille : preision for the approximation of the grid
	# variance to be used in the normal distribution function
	
	## Computes sqrt(n) * sup_x | F_n(x) - Phi(x) | where Phi is the normal distribution function mean 0
	## variance sigma, F_n is the empirical distribution function computed from A

library("mvtnorm")

pval <- function(dis,cv){
	mean(cv>dis)
}

OptParam <- function(X,Y,poids,indice){
	## trouve les param optimaux pour indice = 1: modèle linéaire, 2 : modèle quadratique
	W <- diag(1/abs(X)^poids)
  xf <- matrix(rep(X,each = indice +1),ncol=indice+1,byrow=T)
	ex <- matrix(rep(0:(indice),length(X)),ncol=indice +1 ,byrow=T)
	Xp = xf^ex
	matInv <- solve(t(Xp)%*%W%*%Xp)
	param.optimaux <- rev(matInv%*%t(Xp)%*%W%*%Y)
	return(param.optimaux)	
}
predic <- function(param,x)
{
	xf <- matrix(rep(x,each = length(param)),ncol=length(param),byrow=T)
	ex <- matrix(rep(0:(length(param)-1),length(x)),ncol=length(param),byrow=T)
	par <- matrix(rep(rev(param),length(x)),ncol=length(param),byrow=T)
	Temp <- par * xf^ex
	return(apply(Temp,1,sum))
}

CV <- function(don,indice,poids)
{
	X = rep(don[,1],ncol(don[,-1]))
	Y = c(don[,-1])
	CVc = 0
	corr <- (1/X^poids)
	for(i in 1:length(Y))
	{
		p <- OptParam(X[-i],Y[-i],poids,indice)
		pr <- (predic(p,X[i]) - Y[i])
		CVc = CVc + (pr)^2 * corr[i]
	}
	return(CVc/sum(corr))
}

fit <- function(A,poids,b,nbgrille,stat,indice){
	## A : matrix of data, column 1 contains the covariates
	## poids : 1/X^poids i.e poids = 0 : no additional weight, poids =2 : inverted quadratique weight
	## b : number of bootstrap replicates to use, suggestion is 1000
	## nbgrille : number of points used to approximate the statistic, suggestion is 50
	## stat : 1 is Kolmogorov smirnoff statistic, 2 is cramer von mises
	## indice : 1 is linear, 2 is quadratic 3 is cubic
	
	don <- as.matrix(A)
	nbr <- ncol(A)-1

	W <- diag(rep(1/abs(don[,1])^poids,nbr))	# Weights

# Optimal parameters
	X = rep(don[,1],ncol(don[,-1]))
	Y = c(don[,-1])
	W <- diag(1/abs(X)^poids)
	xf <- matrix(rep(X,each = indice +1),ncol=indice+1,byrow=T)
	ex <- matrix(rep(0:(indice),length(X)),ncol=indice +1 ,byrow=T)
	Xp = xf^ex
	meanXp <- apply(Xp,2,mean)
	matInv <- solve(t(Xp)%*%W%*%Xp)
	matt <- matInv%*%t(Xp)%*%W
	param.optimaux <- rev(matt%*%Y)
	
	#Residuals
	poidsob <- 1/abs(X)^poids
	poidsob <- poidsob/sum(poidsob)
	epsilonp <- (predic(param.optimaux,X) - Y) * sqrt(poidsob)
	var.param <- matInv * var(epsilonp) ## variance of estimated residuals

	## Bootstrap
	n <- length(epsilonp)
	sigma <- sqrt(var(epsilonp))
	
	## Computationi of the data process
	grille = (1:(nbgrille))/(nbgrille+1) ## Approximation over the grid
	grrep <- t(matrix(rep(grille,each=n),byrow=T,ncol=n))
	epsilonrep <- matrix(rep(epsilonp,each=nbgrille),ncol=nbgrille,byrow=T)/sigma
	Fn <- apply(pnorm(epsilonrep)<=grrep,2,mean)
	monproc <- sqrt(n)*(Fn - grille)
	
	proc=matrix(rep(0,b*nbgrille),ncol=b)
	Xb = predic(param.optimaux,X);
	for(i in 1:b)
	{
		ech_boot = sample(epsilonp,n,replace=T)
		Yb = Xb + ech_boot / sqrt(poidsob)
		param.optimauxb <- rev(matt%*%Yb)
		epsilonpb <- (predic(param.optimauxb,X) - Yb) * sqrt(poidsob)
		epsilonrepb <- matrix(rep(epsilonpb,each=nbgrille),ncol=nbgrille,byrow=T)/sqrt(var(epsilonpb))
		Fnb <- apply(pnorm(epsilonrepb)<=grrep,2,mean)
		monprocb <- sqrt(n)*(Fnb - Fn)
		proc[,i] = monprocb
	}

	## Computation of the pvalue
	distr <- rep(0,b)
	mastat <- 0
	if(stat==1){	distr <- apply(abs(proc),2,max)
		mastat <- max(abs(monproc)) }
	if(stat==2){ distr <- apply(proc^2,2,mean)
	 	mastat<- mean(monproc^2) }
	quantileb = quantile(distr,0.95)
	pvalN = mean(distr>mastat)
	
liste <- list(param=param.optimaux,varparam=var.param,pvalN=pvalN,statN=mastat,qBoot=quantileb,esp=epsilonp)
return(liste)
}

fitPlus <- function(A,poids,b,nbgrille,stat)
{
	## Performs the fit plus the partial F text
	fitLin <- fit(A,poids,b,nbgrille,stat,1)
	fitQ <- fit(A,poids,b,nbgrille,stat,2)
	
	SSREGL <- sum(fitLin$esp^2)
	SSREGQ <- sum(fitQ$esp^2)
	degf <- length(fitQ$esp)-3
	stat <- (SSREGL-SSREGQ)/(SSREGQ/(degf))
	pvalF <- 1-pf(stat,1,degf)
	
	L2 <- list(Linear= fitLin,Quadratic = fitQ,Ftest=stat,pvalFtest=pvalF)
	return(L2)
}


