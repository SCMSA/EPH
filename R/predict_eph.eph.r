predict_eph.eph <- function(object, newdata){

  # newdata  : design.grid (e.g coordonnées x,y de chaque point à estimer)
  # coordinates of measurement points
  mesures <-object@x
  # values of measures (observations)
  sorties<-object@y
  # number of input parameters
  NofParam<-object@NofParam
  # number of measures
  NofMeas<-object@NofMeas
  # number of output variables
  NofInk2<-object@NofInk2
  # bounds on intput variable
  boundaries<-object@boundariesin
  # bounds on output variable
  boundariesP<-object@boundariesout
  pas<-object@step
  # number of points to estimate
  NofInk <- nrow(newdata)
  
  save_newdata <- newdata
  NofInk3 <- nrow(save_newdata)
  tempon <- c()
  
  # remove the observations from the points to
  # estimate
  for (i in 1:NofMeas){
    for (j in 1:NofInk){
      compt=0
      for (k in 1:NofParam){
        if ((mesures[i,k]== newdata[j,k]) == TRUE){
          compt=compt+1
          if (compt==NofParam){
            tempon <- c(tempon,j)
          }
        }
      }
    }
  }
  
  temp <- length(tempon) 
  if (temp!= 0){
    for ( i in 1:temp){
      era=temp+1-i
      # in one dimensional case put
      # newdata<-matrix(newdata[-tempon[era],])
      newdata<-newdata[-tempon[era],]
    }
  }
  
  NofParam <- ncol(mesures)
  NofMeas <- nrow(mesures)
  NofInk <- nrow(newdata)

  maxData=rep (0,NofParam)
  minData=rep (0,NofParam)
  x1=matrix (rep(0, NofInk*NofParam), NofInk, NofParam)
  x=matrix (rep(0, NofInk*NofParam), NofInk, NofParam)
  dmax=rep(0,NofParam)
  a=matrix (rep(0, NofMeas*NofParam), NofMeas, NofParam)

  # check if there are NaN in the location of measures
  for ( n in 1:NofMeas) { 
    for (k in 1:NofParam) {
      if (is.na(mesures[n,k])==TRUE){
        stop("sorry, but measure is not full")
      }
      else { 
        if (n == 1) {
          maxData[k] <- mesures[n, k]
          minData[k] <- mesures[n, k]
        }
      } 
    }
  }
  # check NaNs in the value of measures
  for (n in 1:NofMeas){
    for (i in 1:NofInk2){
      if (is.na(sorties[n,i])==TRUE) {
        stop("sorry, but T is not full")
      }
    }
  }
  
  # check NaN in boundaries
  for (j in 1:NofInk){
    for (k in 1:NofParam) {
      if (is.na(boundaries[1, k])==TRUE) {
        stop("sorry, but unknown point is not full")
      }
      else {   
        x1[j,k] <- newdata[j,k]
      }
    }
  }
  
  # compute min and max of input measures coordinates
  for (k in 1:NofParam){
    for (n in 2:NofMeas) {
      if (mesures[n, k] > maxData[k]) {
        maxData[k] <- mesures[n, k]
      }
      if (mesures[n, k] < minData[k]) {
        minData[k] <- mesures[n, k]
      }
    }
  } 
  
  # check if the measures are in the boundaries
  mini <- 1:NofParam
  maxi <- 1:NofParam
  for (j in 1:NofInk){
    for (k in 1:NofParam){
      
      if (boundaries[k+1, 2] < maxData[k]) {
        
        stop( "Sorry, but max for at least one parameter which you chose boundary is smaller than max among measures")
        
      }else{ maxi[k] <- boundaries[k+1 , 2]
      }
      
      if (boundaries[k+1, 1] > minData[k]) {
        stop( "Sorry, but max for at least one parameter which you chose boundary is bigger than min among measures")
        
      }else{ mini[k] = boundaries[k+1 , 1]
      }
    }
  }
  
  tmax<-tmin<-0 
  
  result50=matrix (rep(0, NofInk2*NofInk), NofInk2, NofInk)
  result5=matrix (rep(0, NofInk2*NofInk), NofInk2, NofInk)
  result95=matrix (rep(0, NofInk2*NofInk), NofInk2, NofInk)
  resultprob=matrix (rep(0, NofInk2*NofInk), NofInk2, NofInk)
  esperance=matrix (rep(0, NofInk2*NofInk), NofInk2, NofInk)
  
  save_result50=matrix (rep(0, NofInk2*NofInk3), NofInk2, NofInk3)
  save_result5=matrix (rep(0, NofInk2*NofInk3), NofInk2, NofInk3)
  save_result95=matrix (rep(0, NofInk2*NofInk3), NofInk2, NofInk3)
  save_resultprob=matrix (rep(0, NofInk2*NofInk3), NofInk2, NofInk3)
  save_esperance=matrix (rep(0, NofInk2*NofInk3), NofInk2, NofInk3)
  
  # for each point to estimate
  # for each dimension
  for (m in 1:NofInk2){
    
    tmax = boundariesP[m , 2]
    tmin = boundariesP[m , 1]
    
    vectj<-seq(tmin,tmax,by=pas[m])
    lvectj<-length(vectj)
    
    NofParamValues <- 1+round(round(((tmax - tmin) / pas[m]), 5),0)
    
    proba_density=array (rep(0, NofInk2*NofInk*NofParamValues), c(NofInk2, NofInk, NofParamValues))
    save_proba_density=array (rep(0, NofInk2*NofInk3*NofParamValues), c(NofInk2, NofInk3, NofParamValues))
    
    # for each point to estimate
    for (i in 1:NofInk) {
      
      marqueur<-0
      
      tmax = boundariesP[m , 2]
      tmin = boundariesP[m , 1]
      
      # normalize coordinate of measurement points
      for (k in 1:NofParam){
        for (n in 1:NofMeas) {
          a[n, k] = (mesures[n, k] - mini[k]) / (maxi[k] - mini[k])
        }
        x[i,k] = (x1[i,k] - mini[k]) / (maxi[k] - mini[k])
      }
      
      d=rep(0,NofMeas)
      c1=rep(0,NofMeas)
      lambda<-0
      u1 <- 0
      u <- rep(0,NofMeas)
      
      # compute lambda: propagation parameter
      for (n in 1:NofMeas){
        for (k in 1:NofParam){
          if (abs(1 - a[n, k]) >= abs(a[n, k])) {
            dmax[k] <- (1 - a[n, k]) ^ 2
          }else{ dmax[k] <- (a[n, k]) ^ 2
          }
          u[n] <- u[n] + dmax[k]
        }
        nju <- 1 + round(round(((tmax - tmin) / pas[m]), 5),0)
      }
      u<-sqrt(u)
      lambda <- log(nju)[1] / (max(u))
      
      # compute distance between each measures and point to estimate
      for (n in 1:NofMeas){
        for (k in 1:NofParam){
          d[n] <- (x[i,k] - a[n, k]) ^ 2 + d[n]
        }
      }
      
      # calcul de sum_{i} d_{i}
      dsum <- 0
      for (n in 1:NofMeas){
        if (identical(d[n], 0) == TRUE){
          #    print ('The unknown point has the same parameters as ' , n , '-th measure point')
          Expect <- sorties[n,m]
          
        }else{ dsum <- dsum + d[n] ^ (-NofParam / 2)
        }
      }
      
      # compute normalisation constant
      c1=rep(0,NofMeas)
      for (n in 1:NofMeas){
        for (j in seq(tmin,tmax,by=pas[m])){
          j <- round(j, 5)
          c1[n] <- exp(-((pi * ((j - sorties[n,m]) ^ 2)) / (pas[m] ^ 2)) * exp(1 - 2 * lambda * (d[n] ^ 0.5)) + 0.5 - lambda * (d[n] ^ 0.5)) + c1[n]
        }
      }
      
      # compute probability law
      # the probability law generated by each observation n, p_{n,j}(X)
      # are combinated : the weight to combinate depends on the distance to the point
      # to reconstruct X      
      p=rep(0,lvectj)
      l <- 1
      sum1 <-0
      for (j in seq(tmin,tmax,by=pas[m])){
        j <- round(j, 5)
        for (n in 1:NofMeas){
          p[l] <- (1 / c1[n]) * (exp(-((pi * ((j - sorties[n,m]) ^ 2)) / (pas[m] ^ 2)) * exp(1 - 2 * lambda * (d[n] ^ 0.5)) + 0.5 - lambda * (d[n] ^ 0.5)) * (1 / (d[n] ^ (NofParam / 2))))/dsum + p[l]
        }
        sum1 <- p[l] + sum1
        l <- l + 1
      }
      # save probability law
      save_proba_density[m,i,]=p  
      
      # compute expectation
      Expect <- 0
      h <- 1
      for (j in seq(tmin,tmax,by=pas[m])){
        j <-  round(j, 5)
        Expect <-  (j * (p[h])) + Expect
        h <-  h + 1
      }
      esperance[m,i]=Expect
      
      # compute mediane
      stock<-rep(0, 2)
      stock[1]<- tmin
      stock[2]<- 0
      h<-1
      for (j in seq(tmin,tmax,by=pas[m])){
        j <-  round(j, 5)
        if (stock[2]<0.5){
          stock[1] = j
          stock[2] = stock[2] + (p[h])
        } 
        h <-  h + 1 
      }
      result50[m,i]=stock[1]
      
      # compute most probable value
      stock<-rep(0, 2)
      stock[1]<- tmin
      stock[2]<- 0
      h<-1
      for (j in seq(tmin,tmax,by=pas[m])){
        j <-  round(j, 5)
        if (stock[2]<p[h]){
          stock[1] = j
          stock[2] = (p[h])
        }
        h <-  h + 1 
      }
      resultprob[m,i]=stock[1]
      
      # compute quantile 5
      stock<-rep(0, 2)
      stock[1]<- tmin
      stock[2]<- 0
      h<-1
      for (j in seq(tmin,tmax,by=pas[m])){
        
        j <-  round(j, 5)
        if (stock[2]<0.025){
          stock[1] = j
          stock[2] = stock[2] + (p[h])
        }
        h <-  h + 1 
      }
      result5[m,i]=stock[1]
      
      # compute quantile 95
      stock<-rep(0, 2)
      stock[1]<- tmin
      stock[2]<- 0
      h<-1
      for (j in seq(tmin,tmax,by=pas[m])){
        
        j <-  round(j, 5)
        if (stock[2]<0.975){
          stock[1] = j
          stock[2] = stock[2] + (p[h])
        }
        h <-  h + 1 
      }
      result95[m,i]=stock[1]
    }
  }

  y=cbind(y)
  NofMeas <- nrow(mesures)
  
  # compute features for each points to reconstruct
  # including the observation points
  # if they are to reconstruct
  save_proba_density_temp=save_proba_density
  for (m in 1:NofInk2){
    h <- 1
    for (i in 1:NofInk3){
      Expect <- 0 
      ok <- 0
      
      # if the point to reconstruct is on
      # an observation point
      for (g in 1:dim(mesures)[1]){
        if (all(save_newdata[i,]==mesures[g,])){
          ok <- 1 
          save_result5[m,i]=sorties[g,]
          save_result50[m,i]=sorties[g,]
          save_result95[m,i]=sorties[g,]
          save_resultprob[m,i]=sorties[g,]
          # the probability law is a Dirac centered in the 
          # value of the observation
          save_proba_density[m,i,]=matrix(rep(0,nju[1]),nju[1])
          save_proba_density[m,i,round(sorties[g,]/pas[m]-tmin)] = 1
          save_esperance[m,i]=sorties[g,m]
        }
      }
      # otherwise : for all other points to reconstruct
      # save all features
      if (ok==0){ 
        save_result5[m,i]=result5[m,h]
        save_result50[m,i]=result50[m,h]
        save_result95[m,i]=result95[m,h]
        save_resultprob[m,i]=resultprob[m,h]
        save_esperance[m,i]=esperance[m,h]
        save_proba_density[m,i,]=save_proba_density_temp[m,h,]
        h <- h + 1
      } 
    }
  }
  
  output.list <- list()
  output.list$lower95 <- save_result5
  output.list$mean <- save_result50
  output.list$upper95 <- save_result95
  output.list$probable <- save_resultprob
  output.list$esperance <- save_esperance 
  output.list$proba_density <- save_proba_density 
  return(output.list)
  
  
}
