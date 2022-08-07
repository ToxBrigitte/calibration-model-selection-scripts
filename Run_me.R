####################
## File to be run ##
####################

## Parameters to be set
directory <- "~/Desktop/Bridge/Calibration" # set work directory : folder where the calibration data is
																						# and the two R sheets Calibration2Ks and Calibration2CVM
filename <- "Model.txt" # name of the file that contains the raw data (calibration data)
result_filname <- "My_Results.txt" # name of the file that contains the results
dec <- "." # , if the decimals are ex. 0,5 . if the decimals are ex 0.5
nbgrille = 50 # global paramter for integral approximations
b <- 1000 #global paramter for bootstraping
stat <- 2 # 1: KS, 2: CVM
alpha <- 0.05 ## The alpha to use in hypothesis testing
## end of parameter setting


#################################################################################################
####################################### Code ####################################################
#################################################################################################
source("CodeV5.R") ## Contains the main coded functions
A <- read.table(filename,header=F,dec=dec) ## reading the data
A <- as.matrix(A) 
stats = c("Kolmogorov Smirnov","Cramer von Mises")


Decision <- rep(0,3)
DCV <- rep(0,3)

# weight : 0
weight <- 0;
P0 = fitPlus(A,weight,b,nbgrille,stat)
coP0 <- matrix(c(rev(P0$Linear$param),0,rev(P0$Quadratic$param)),ncol=3,byrow=T) ## Contains the fitted coefficients
colnames(coP0) <- c("b0","b1","b2")
rownames(coP0) <- c("Linear","Quadratic")


List0 <- list(Fitted_Param = coP0, Pval_Lin = P0$Linear$pvalN,Stat_Lin = P0$Linear$statN ,
			Pval_Quad = P0$Quadratic$pvalN,Stat_Quad = P0$Quadratic$statN, Pval_PartialF = P0$pvalFtest)

Decision[1] = 1 + (List0$Pval_PartialF<0.05)

# weight : 1
weight <- 1;
P1 = fitPlus(A,weight,b,nbgrille,stat)
coP1 <- matrix(c(rev(P1$Linear$param),0,rev(P1$Quadratic$param)),ncol=3,byrow=T) ## Contains the fitted coefficients
colnames(coP1) <- c("b0","b1","b2")
rownames(coP1) <- c("Linear","Quadratic")


List1 <- list(Fitted_Param = coP1, Pval_Lin = P1$Linear$pvalN,Stat_Lin = P1$Linear$statN ,
	Pval_Quad = P1$Quadratic$pvalN,Stat_Quad = P1$Quadratic$statN, Pval_PartialF = P1$pvalFtest)

Decision[2] = 1 + (List1$Pval_PartialF<0.05)

# weight : 2
weight <- 2;
P2 = fitPlus(A,weight,b,nbgrille,stat)
coP2 <- matrix(c(rev(P2$Linear$param),0,rev(P2$Quadratic$param)),ncol=3,byrow=T) ## Contains the fitted coefficients
colnames(coP2) <- c("b0","b1","b2")
rownames(coP2) <- c("Linear","Quadratic")


List2 <- list(Fitted_Param = coP2, Pval_Lin = P2$Linear$pvalN,Stat_Lin = P2$Linear$statN ,
	Pval_Quad = P2$Quadratic$pvalN,Stat_Quad = P2$Quadratic$statN, Pval_PartialF = P2$pvalFtest)

Decision[3] = 1+ (List2$Pval_PartialF<0.05)

## Variance Test for weight selection
Ano  = A[,-1]/sqrt(A[,1])
spoids = sum(1/sqrt(A[,1]))
var1 = apply(Ano/spoids,1,var)
Ano2 = A[,-1]/A[,1]
spoids2 = sum(1/(A[,1]))
var2 = apply(Ano2/spoids2,1,var)
var0 = apply(A[,-1]/length(A[,1]),1,var)
s1 = signif(sum((var0-mean(var0))^2),3)
s2 = signif(sum((var1-mean(var1))^2),3)
s3 = signif(sum((var2-mean(var2))^2),3)
sss = c(s1,s2,s3)
z1 = paste("Variance test for weight selection" )
z2 = paste("Scores:" , "No weight:",s1,"x^(-1):",s2, " x^(-2):",s3 )
ccc = c("no weight","x^(-1)","x^(-2)")
ccb = paste("Selected weight: ",ccc[which(sss==min(sss))],sep="")
zaa = paste(z1,z2,ccb,sep="\n")

best <- which(sss==min(sss))

awiner = eval(parse(text = paste("List",best-1,sep="")))

## Ftest Hétéro
pop1 = A[1,-1]
pop2 = A[nrow(A),-1]
ttt = var.test(pop1,pop2,alternative = "less")
bb = paste("F-test for heteroscedasticity")
bb1 = paste("p-value: ",signif(ttt$p.value,4))
ddd = "No"
if(ttt$p.value<alpha ){ddd="Yes"}
bb2 = paste("Weighting needed:",ddd,sep=" ")
popl = paste(bb,bb1,bb2,sep="\n")


v = paste("Partial F-test for model order selection")
v1 = paste("p-value",signif(awiner$Pval_PartialF,4),sep=": ")
v2 = paste("Model selected: ",c("linear","quadratic")[Decision[best]],sep="")
vv = paste(v,v1,v2,sep="\n")


w = paste("Normality of the standerdized residuals")
w1 = paste("Test used: ",stats[stat],sep="")
aw= eval(parse(text = paste("awiner$Pval_",c("Lin","Quad")[Decision[best]],sep="")))
w2 = paste("p-value: ",signif(aw,4))
ww = "No"
if( alpha<aw)
{
	ww = "Yes"
}
ww1 = paste("Validation test passed:",ww,sep=" ")
vw =paste(w,w1,w2,ww1,sep="\n")

md = c("Linear","Quadratic")[Decision[best]]
c("1","1/x","1/x^2")[best]
if(ww=="Yes")
{
zw=paste("Model selected: ",md,", ",c("1","1/x","1/x^2")[best],sep="" )
zw2 = paste("Calibration equation:")
if(md=="Linear")
{
	co1 = awiner$Fitted_Param[1,1]
	co2 = awiner$Fitted_Param[1,2]
	zw1 = paste(signif(co2,4)," x + ",signif(co1,4),sep="")
}
if(md=="Quadratic")
{
	co1 = awiner$Fitted_Param[2,1]
	co2 = awiner$Fitted_Param[2,2]
	co3 = awiner$Fitted_Param[2,3]
	zw1 = paste(signif(co3,4)," x^2 + ",signif(co2,4)," x + ",signif(co1,4),sep="")
}
vw2 = paste(zw,zw2,zw1,sep="\n")
}
if(ww == "No"){
	vw2 = paste("No model selected, validation test failed")
}

## Plots
namePlot = strsplit(result_filname,".",fixed=TRUE)[[1]][1]
# Plot of variance
var_level = apply(A[,-1],1,var)
pdf(file=paste(namePlot,"Variance.pdf",sep="_"))
plot(A[,1],var_level,xlab="Concentration",ylab="Variance",main="Variance plot",mgp = c(2, 0.8, 0), axes = T)
dev.off()
## Calibration curve
f <- function(x)
{
	predic(rev(awiner$Fitted_Param[Decision[best],]),x)
}
Cal.dots = predic(rev(awiner$Fitted_Param[Decision[best],]),A[,1])
pdf(file=paste(namePlot,"Calibration curve.pdf",sep="_"))
plot(rep(A[,1],each=(ncol(A[,-1]))),c(t(A[,-1])),xlab="Concentration",ylab="Signal",main="Calibration curve",mgp = c(2, 0.8, 0), axes = T)
curve(f,add=T)
dev.off()
mt = paste(popl,zaa,vv,vw,vw2,sep="\n \n")
write(mt,file=result_filname)
cat(mt)

