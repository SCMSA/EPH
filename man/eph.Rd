% File man/eph.Rd
\name{eph}
\alias{eph}
\title{builds eph object}
\description{
builds a EPH object
}
\arguments{
	\item{Design}{Points to reconstruct:
			design.grid is the location of the points to reconstruct, 
			design.fact is the location of the point}
	\item{Response}{Observation point: 
			response.grid is the expected (real) value to the points to reconstruct, 
			response.fact is the value of each observations}
	\item{Boundaries}{value of the measure at each observation point}
	\item{BoundariesP}{value of the measure at each observation point}
	\item{Step}{value of the measure at each observation point}
}
\value{
	\item{x}{coordinates of each observation point}
	\item{y}{value of the measure at each observation point}
	\item{NofParam}{number of input parameter}
	\item{NofMeas}{number of observation points}
	\item{NofInk2}{number of output variables}
	\item{Boundariesin}{value of min and max for each input parameter}
	\item{Boundariesout}{value of min and max for each output variables}
	\item{Step}{discretisation step for each input and output parameter}
}
\author{Gottfried Berton}
