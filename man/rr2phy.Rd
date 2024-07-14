\name{rr2phy}
\alias{rr2phy}
\docType{package}
\title{
Computes Phylogenetic eigen Vectors Regression and evalute the effects of evolution and environment on plant traits using lambda method.
}
\description{
the varation decoposition of phylogeny and adapation enviroment condition
}
\details{
\tabular{ll}{
Package: \tab rr2PVR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2023-06-07\cr
License: \tab GPL >= 2\cr
}

}
\author{
Wuchao gao; Qi li
  
Maintainer: Wuchao gao <sichuangwc@gmail.com>
}
\references{
Diniz-Filho, J.A.F., Sant'Ana, C.E.R. and Bini, L.M. (1998). An eigenvector method for estimating phylogenetic inertia. Evolution, 52:1247-1262.

Keck, F., Rimet, F., Bouchez, A. and Franc, A. (2016), phylosignal: an R package to measure, test, and explore the phylogenetic signal. Ecol Evol, 6: 2774-2780. https://doi.org/10.1002/ece3.2051
}

\keyword{ package }
\seealso{

\code{\link{phyPVR}}
}
\examples{
library(rr2phy)  
library(rdacca.hp)  
library(ape)   
library(PVR)   
library(phylosignal)   
library(phylobase)   
##Evaluate the relative importance of species trait in phylogeny and a single environment factor  
tree <- rcoal(10)  
x <- PVRdecomp(tree)  
trait <-  rnorm(10,1,0.01)  
envvar <- rnorm(10,2,0.05)  
phyPVR(x, trait = trait, envVar = envvar, method = "Moran.I")    
#test phylogenic siginal     
exGeo1<-phylo4d(tree,tip.data =trait )    
phyloSignal(exGeo1)    
##Evaluate the relative importance of species trait in phylogeny and multiple environment factors  
tree <- rcoal(10)  
x <- PVRdecomp(tree) 
trait <-  rnorm(10,1,0.01)  
envvar <-  matrix(rnorm(100,2,0.1), nrow = 10, ncol = 10)  
phyPVR(x, trait = trait, envVar = envvar, method = "Moran.I")  
#test phylogenic siginal     
exGeo1<-phylo4d(tree,tip.data =trait )    
phyloSignal(exGeo1)    
##test real data   
library(caper)   
data(shorebird)   
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)   
x <- PVRdecomp(shorebird$phy)   
envvar <- matrix(rnorm(71,2,0.1), nrow = 71, ncol = 10)   
phyPVR(x, trait = shorebird$data$M.Mass, envVar = envvar, method = "Moran.I")   
#test phylogenic siginal        
exGeo1<-phylo4d(shorebird$phy,tip.data = shorebird$data$M.Mass)   
phyloSignal(exGeo1)   
#barplot.phylo4d(exGeo1)   
}
