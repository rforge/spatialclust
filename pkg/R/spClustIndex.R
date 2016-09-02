#' Cluster Validity Index
#'
#' @description This function used to validate the clustering result
#'
#' @param fgwc result(object) from fgwc clustering
#'
#' @return validity indeks
#'
#' @examples
#' #load data example
#' X <- example
#'
#' #if using matrix distance
#' distance <- dist
#'
#' #if using shapefile
#' #library(rgdal) for call readOGR
#' #distance <- readOGR(dsn = 'folder/.',"shapefile name")
#'
#' #load population data
#' pop <- population
#'
#' clust <- fgwc(X,pop,distance,K=2,m=1.5,beta=0.5)
#'
#' #show cluster validation
#' spClustIndex(clust)
#'
#'
#' @seealso \code{\link{visualize}} for cluster visualizatiion
#' \code{\link{scale}} for data scalling
#'
#' @export

spClustIndex <- function(fgwc){
  partitionMatrix <- fgwc$U
  matriksPusatKlaster <- fgwc$V
  data <- fgwc$data
  m <- fgwc$m

  C <- ncol(partitionMatrix)
  N <- nrow(partitionMatrix)
  D <- ncol(data)
  dist <- matrix(0,N,C)
  distV <- matrix(0,C,C)

  #hitung jarak
  for(i in 1:N){
    for(k in 1:C){
      pembilang <- 0
      for(j in 1:D){
        pembilang <- pembilang + abs(data[i,j] - matriksPusatKlaster[k,j]) ^ 2
      }
      dist[i,k] <- pembilang
    }
  }

  #jarak antar klaster
  for(k1 in 1:C){
    for(k2 in 1:C){
      distV[k1,k2]<-t(matriksPusatKlaster[k1,]-matriksPusatKlaster[k2,])%*%(matriksPusatKlaster[k1,]-matriksPusatKlaster[k2,])
    }
  }

  PC <- 0
  for(i in 1:C){
    for(j in 1:N){
      PC <- PC + (partitionMatrix[j,i]^2)
    }
  }
  PC <- PC/N

  CE <- 0
  for(i in 1:C){
    for(j in 1:N){
      CE <- CE + partitionMatrix[j,i]*log(partitionMatrix[j,i],base = exp(1))
    }
  }
  CE <- -1*CE/N

  SC <- 0
  pembilang<-matrix(0,N,C)
  for(i in 1:N)
  {
    for(k in 1:C){
      pembilang[i,k]<-(t(data[i,]-matriksPusatKlaster[k,])%*%(data[i,]-matriksPusatKlaster[k,]))*(partitionMatrix[i,k]^m)
    }
  }

  penyebut <- 0
  SC <- 0
  for(i in 1:C){
    penyebut <- sum(partitionMatrix)*sum(distV)
    SC <- SC+(sum(pembilang[,i])/penyebut)
  }

  S <- 0
  pembilang<-matrix(0,N,C)
  for(i in 1:N)
  {
    for(k in 1:C){
      pembilang[i,k]<-(t(data[i,]-matriksPusatKlaster[k,])%*%(data[i,]-matriksPusatKlaster[k,]))*(partitionMatrix[i,k]^2)
    }
  }

  penyebut<-min(distV[distV>0])*N
  S<-sum(pembilang)/penyebut

  XB <- 0
  pembilang<-matrix(0,N,C)
  for(i in 1:N)
  {
    for(k in 1:C){
      pembilang[i,k]<-(t(data[i,]-matriksPusatKlaster[k,])%*%(data[i,]-matriksPusatKlaster[k,]))*(partitionMatrix[i,k]^m)
    }
  }

  penyebut<-min(distV[distV>0])*N
  XB<-sum(pembilang)/penyebut


  dimensi <- dim(data)
  n <- dimensi[1]
  d <- dimensi[2]
  c <- ncol(partitionMatrix)

  #deviation
  temp1 <- 0
  for(j in 1:c){
    temp <- 0
    for(k in 1: n){
      temp <- temp + abs(data[k,] - matriksPusatKlaster[j,])^2
    }
    temp1 <- temp1 + temp / n
  }
  dev <- sum(temp1) / c


  #sd max
  for(j in 1:c){
    for(k in 1:c){
      if(j != k){
        temp <- temp + abs(matriksPusatKlaster[k,] - matriksPusatKlaster[j,])^2
        jumlah <- sum(temp)
        if(j == 1){
          max <- jumlah
        }else if(max < jumlah){
          max <- jumlah
        }
      }
    }
  }
  sd <- max

  #ifv
  temp2 <- 0
  for(j in 1:c){
    temp1 <- 0
    temp <- 0
    for(k in 1:n){
      temp <- temp + log(partitionMatrix[k,j], base=2)
    }

    for(k in 1:n){
      temp1 <- temp1 + (partitionMatrix[k,j]^ 2 * (log(c, base=2) - (temp/n))^2)
    }
    temp2 <- temp2 + ((temp1 / n)* (sd / dev))
  }

  ifv <- (temp2 / c)

  hasil<-c(PC,CE,S,XB,ifv)
  names(hasil)<-c("PC index","CE index","Separation index","XB index","IFV index")
  return(hasil)

}
