% File man/predict_eph.Rd
\name{predict_eph.eph}
\alias{predict_eph.eph}
\title{run EPH}
\arguments{
	\item{object}{object of class EPH}
	\item{newdata}{coordinates on which performing the estimation}
}
\value{
	A list of vectors, the components are :
	\item{mean}{quantile 50 of the probability law}
	\item{upper95}{quantile 95 of the probability law}
	\item{probable}{most probable value}
	\item{esperance}{expected value of the probabilty law}
	\item{lower95}{quantile 5 of the probability law}
	\item{proba_density}{raw probabilty distribution computed by the EPH}
}
\references{Olga Zeydina et Bernard Beauzamy : Probabilistic Information Transfer. 
Societe de Calcul Mathematique SA. 
ISBN: 978 2 9521458 6 2, ISSN : 1767 1175. Reli, 208 pages, mai 2013.

Reactor Safety and Incomplete Information: Comparison of Extrapolation Methods for the Extension of 
Computational CodeProceedings of ICAPP 2015 Nice, France, May 3-6, 2015
Paper 15377
F. GODAN, O. ZEYDINA, B. BEAUZAMY (Societe de Calcul Mathematique)
Y. RICHET (Institut de Radioprotection et de Surete Nucleaire)
}
\description{
this function estimates the expected value, the mean and mediane value of the 
probability law returned by the EPH in each point of the domain.
At the observation points, the algorithm returns the real value.
}
\author{Gottfried Berton}
\examples{
 
## a 2D example : reconstruct the branin function
library(DiceKriging)

d <- 2; 
#number of observations
n <- 9

# bound output variable (e.g temperature)
tmax<- (500) 
tmin<- (0)

# position of observations
design.fact <- expand.grid(x1=seq(0,1,length=3), x2=seq(0,1,length=3))
# generate value of observations
y <- apply(design.fact, 1, branin)
response.fact<-cbind(y)

setClass("eph",slots=c(NofParam="integer", NofMeas="integer",
NofInk2="integer", x="matrix", y="matrix", step="matrix", 
boundariesin="matrix",boundariesout="matrix", covariance="covTensorProduct"))
 
# bound of the parameter
mat=matrix(data=c(tmin,0,0,tmax,1,1), ncol=2)
boundariesP=rbind(c(tmin,tmax))

pas=cbind(c(0.1,0.1,0.1))
eph1 <- eph(design.fact, response.fact, mat, boundariesP, pas)

# generate grid points to reconstruct
n.grid <- 5
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x1=x.grid, x2=y.grid)
response.grid <- apply(design.grid, 1, branin)
nb_to_estimate<-nrow(design.grid)

# run EPH method
m<-1
eph1 <- eph(design.fact, response.fact, mat, boundariesP, pas)
system.time(carac <- predict_eph.eph(eph1, design.grid))
loi <- carac$proba_density
esp <- carac$esperance
 
}
