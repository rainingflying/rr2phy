# rr2phy
the rr2phy packages mianly provides the Moran.I method for selecting the optimal phylogenetic eigenvectors and allows for the direct evaluation of the effects of evolution and environment on plant traits.The phylogenetic eigenvector regression (PVR) starts by performing an eigendecomposition of a pairwise double-centered phylogenetic distance matrix between species. The eigenvectors estimated values express phylogenetic trends in data and residuals express independent evolution of each species. 
## To install the rr2phy packages from Github using devtools:
library(devtools)     
devtools::install_github("rainingflying/rr2phy")  
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
##Evaluate the relative importance of species trait in phylogeny and multiple environment factors  
tree <- rcoal(10)  
x <- PVRdecomp(tree) 
trait <-  rnorm(10,1,0.01)  
envvar <-  matrix(rnorm(100,2,0.1), nrow = 10, ncol = 10)  
phyPVR(x, trait = trait, envVar = envvar, method = "Moran.I")  
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
