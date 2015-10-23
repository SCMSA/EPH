# Package-EPH
 Authors
--------------------------------------------------------

Gottfried Berton, gottfried.berton at scmsa.eu

Installation
-----------------------------------------------------------------------

You can install the latest version of the code using the devtools R package.

# Install devtools, if you haven't already.
install.packages("devtools")

library(devtools)
install_github("EPH", "SCMSA")


Why should I use it ?
-----------------------------------------------------------------------

The 'EPH' package runs the experimental probabilistic hypersurface method. This method aim at
reconstructing or predicting values based on available data. The principle is the propagation of the
information at observation points toward unknown points. The result obtained for each point to reconstruct
is a probability density. The importance given to the input data is based on their distance to the point
to reconstruct. The more a measure is close to the point to reconstruct, the greater will be its
influence on the probability density.

How do I use it ?
-----------------------------------------------------------------------

The dimension of space is chosen by the user. The number of observation points is arbitrary too.
The user set the following parameters :

- the min and max for each input parameter, and output variable
- value of the measures for each observation point
- the coordinates of the points to reconstruct
- discretisation step for the different parameters (input and output)
 
What does the package includes ?
-----------------------------------------------------------------------

- the function class eph that creates an object of class eph
- the function predict_eph.eph that estimates the expected value of the probability law, as well as the mediane, quantiles, etc.
  for each point to reconstruct.
  

