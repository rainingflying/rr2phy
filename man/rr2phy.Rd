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
library(PVR)
library(ape)
library(rdacca.hp)
library(phylosignal)
library(rr2phy)
tree <- rcoal(10)
x <- PVRdecomp(tree)
trait <-  rnorm(10,1,0.01)
envvar <- rnorm(10,2,0.05)
phyPVR(x, trait = trait, envVar = envvar, method = "lambda")
}
