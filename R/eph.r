eph <- function (design, response, boundaries, boundariesP, step)
{
 
  model<-new("eph")
  data<-data.frame(design)
  x <-as.matrix(data)
  y=as.matrix(response)
  model@x <- x
  model@y<- y
  model@NofParam <- ncol(x)
  model@NofMeas <- nrow(x)
  model@NofInk2 <-ncol(y)
  model@boundariesin <- boundaries
  model@boundariesout <- boundariesP
  model@step <- step
  known.covparam <- "None"
  covtype="matern5_2"
 
  # loi de probabilité des incertitudes des valeures de chaque mesures n
  # ptheta est une matrice de taille (nombre de mesures X nombre d'intervalles de discretisation sur le param de sortie)
  
  return(model)
}
