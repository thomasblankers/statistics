## this function was developped with Dean C. Adams, Professor at Iowa State.
## the purpose was to perform multivariate phylogenetic least squares analysis on size corrected morphometric measurements


#Burnaby's size-correction
#Find isometric size vector
f.iso = array(1,dim=ncol(data))/sqrt(ncol(data))
# calculate shape based on isometric size vector
shape = log(data)%*%(diag(ncol(data))-(f.iso%*%solve(t(f.iso)%*%f.iso)%*%t(f.iso)))
#perform pca
shape.pca = prcopm(shape, scale.=TRUE)

# Multivariate pGLS
#For Brownian motion model

library(ape)
#phy is a phylogeny in ape's phylo format,
#y.mat is a matrix with response variables (e.g. shape morphology),
#x.mat is a matrix with predictor variables (e.g. microhabitat measurements)
mult.pgls<-function(phy,y.mat,x.mat){
phy.mat<-vcv(phy)
x.mat<-cbind(matrix(1,length(phy$tip.label)),x.mat)
x.mat<-x.mat[rownames(phy.mat),]
y.mat<-y.mat[rownames(phy.mat),]
mDnew<-solve(svd(phy.mat)$u%*%diag(sqrt(svd(phy.mat)$d))%*%t(svd(phy.mat)$u))
ynew<-mDnew %*% y.mat
xnew<-mDnew %*% x.mat
summary=summary(manova(lm(ynew~xnew-1)))
AIC=extractAIC(lm(ynew~xnew-1), k=2, scale=0)
return(list(summary=summary, AIC=AIC))
}

# For Ornstein-Uhlenbeck model

library(ape)
#phy is a phylogeny in ape's phylo format,
#y.mat is a matrix with response variables (e.g. shape morphology),
#x.mat is a matrix with predictor variables (e.g. microhabitat measurements)
library(geiger)
alpha=fitContinuous(phy, x.mat, method=”OU”)$alpha
phy.ou=corMartins(alpha, phy)

mult.pgls<-function(phy.ou,y.mat,x.mat,phy) {
phy.mat<-vcv(phy.ou)
x.mat<-cbind(matrix(1,length(phy$tip.label)),x.mat)
x.mat<-x.mat[rownames(phy.mat),]
y.mat<-y.mat[rownames(phy.mat),]
mDnew<-solve(svd(phy.mat)$u%*%diag(sqrt(svd(phy.mat)$d))%*%t(svd(phy.mat)$u))
ynew<-mDnew %*% y.mat
xnew<-mDnew %*% x.mat
summary=summary(manova(lm(ynew~xnew-1)))
AIC=extractAIC(lm(ynew~xnew-1), k=2, scale=0)
return(list(summary=summary, AIC=AIC))
}
