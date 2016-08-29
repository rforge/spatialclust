#'  Fuzzy Geographically Weighted Clustering (FGWC) optimized by Gravitational Search Algorithm
#'
#' @description This function used to perform Fuzzy Geographically Weighted Clustering of X dataset.
#' by using this function the initialization phase of FGWC will be optimized using Gravitational Search Algorithm
#'
#' @param X data frame n x p
#' @param population dataset 1 x n number of population each region (row)
#' @param distance shapefile or distance matrik n x n
#' @param K specific number of cluster (must be >1)
#' @param m fuzzifier / degree of fuzziness
#' @param beta proportion of geographically effect (if 0 equal Fuzzy C-Means)
#' @param a power for increase population effect
#' @param b power for increase distance effect
#' @param max.iteration maximum iteration to convergence
#' @param threshold threshold of convergence
#' @param RandomNumber specific seed
#'
#' @return func.obj objective function that calculated.
#' @return U matrix n x K consist fuzzy membership matrix
#' @return V matrix K x p consist fuzzy centroid
#' @return D matrix n x K consist distance of data to centroid that calculated
#' @return Clust.desc cluster description (dataset with additional column of cluster label)
#'
#' @details This function perform Fuzzy Geographically Weighted Clustering by G.A Mason and R.Jacobson (2007).
#' Fuzzy Geographically Weighted Clustering is one of fuzzy clustering methods to clustering dataset
#' become K cluster. Number of cluster (K) must be greater than 1. To control the overlaping
#' or fuzziness of clustering, parameter m must be specified.
#' Maximum iteration and threshold is specific number for convergencing the cluster.
#' Random Number is number that will be used for seeding to firstly generate fuzzy membership matrix.
#' population dataset, shapefile or distance matrix is used to give geographically weighted for membership matrix.
#' @details Clustering will produce fuzzy membership matrix (U) and fuzzy cluster centroid (V).
#' The greatest value of membership on data point will determine cluster label.
#' Centroid or cluster center can be use to interpret the cluster. Both membership and centroid produced by
#' calculating mathematical distance. Fuzzy Geographically Weighted Clustering calculate distance with Euclideans norm. So it can be said that cluster
#' will have sperichal shape of geometry.
#'
#' @references G. A. Mason and R. D. Jacobson.(2007). Fuzzy Geographically Weighted Clustering, in Proceedings of the 9th International Conference on Geocomputation, no. 1998, pp. 1-7
#' @references Bezdek, J. C., Ehrlich, R., & Full, W. (1984). FCM: The Fuzzy C-Means Clustering Algorithm. Computers and Geosciences Vol 10, 191-203
#'
#' @export
#'

fgwc.gsa<- function(X,population,distance,K=2,m=2,beta=0.5,a=1,b=1,max.iteration=100,threshold=10^-5,
                RandomNumber=0) {
  ## Set data
  data.X <- as.matrix(X)
  population <- as.matrix(population)
  n <- nrow(data.X)
  p <- ncol(data.X)
  alfa <- 1- beta
  map <- NULL

  if(is.matrix(distance)){
    distance <- as.matrix(distance)
  }else{
    library(rgeos)
    centroid <- gCentroid(distance,byid = T)
    map <- distance
    distance <- as.matrix(spDists(centroid, longlat=T))
  }

  ##Initiation Parameter##
  if ((K <= 1) || !(is.numeric(K)) || (K %% ceiling(K) > 0))
    K = 2
  if ( (m <= 1) || !(is.numeric(m)))
    m = 2
  if (RandomNumber > 0)
    set.seed(RandomNumber)

  ## Membership Matrix U (n x K)
  U <- matrix(runif(n * K,0,1),n,K)

  ## Prerequirement of U:
  ## Sum of membership on datum is 1
  U <- U / rowSums(U)

  ## Centroid Matrix V (K x p)
  V <- matrix(0,K,p)

  ## Distance Matrix
  D <- matrix(0,n,K)

  #GSA
  pop <- nrow(data.X)
  const <- 650

  #velocity
  Velocity <- matrix(0,n,K)
  Vm <- Velocity
  force <- matrix(0,K,pop)

  #candidate solution
  fitness <- vector("numeric")
  CurrentFitness <- vector("numeric")
  best <- vector("numeric")
  worst <- vector("numeric")
  mass <- vector("numeric")
  Mass <- vector("numeric")

  arrayBestfit <- vector("numeric")
  fit <- vector("numeric")


  #weighted membership matrix
  mi.mj <- (population %*% t(population)) ^ b
  diag(distance) <- Inf
  weighted <- mi.mj / distance

  #give weighted for membership
  membership <- weighted%*%U
  summ <- as.matrix(rowSums(membership))
  #bagi tiap baris matriks U dengan pembagi
  for(j in 1:K){
    for(i in 1:n){
      U[i,j] <- (alfa*U[i,j])+(beta*(membership[i,j]/summ[i,1]))
    }
  }

  noAgen <- 0
  t <- 0
  while(t < max.iteration){
    t <- t + 1
    #acceleration
    A <- matrix(0,nrow = pop,ncol = K)

    #best fitness (minimal)
    fmin <- Inf
    for(iter in 1:pop){
      #jumlah agen
      noAgen <- (t-1)*pop+iter
      G <- 1*exp(-5*iter)

      #matrik partisi
      U.old <- U

      #calculate centroid
      V <- t(U.old ^ m) %*% data.X / colSums(U.old ^ m)
      for (k in 1:K)
      {
        #Distance calculation
        for (i in 1:n)
        {
          D[i,k] = t(data.X[i,] - V[k,]) %*%
            (data.X[i,] -V[k,])
        }
      }

      #hitung jarak data ke pusat klaster
      f <- sum(U.old ^ m * D)
      fit[iter] <- f
      CurrentFitness[noAgen] <- fit[iter]
      fitness[noAgen] <- fit[iter]

      if(fmin > fitness[noAgen]){
        fmin <- fitness[noAgen]
      }

      arrayBestfit[noAgen] <- fmin

      if(noAgen > 1){
        if(abs(arrayBestfit[noAgen] - arrayBestfit[noAgen-1]) < threshold){
          break
        }
      }

      best[noAgen] <- min(CurrentFitness)
      worst[noAgen] <- max(CurrentFitness)

      for(i in 1:noAgen){
        if( i == 1){
          mass[noAgen] <- 0
        }else{
          mass[i] <- (CurrentFitness[i]  -   worst[noAgen])/(best[noAgen]-worst[noAgen])
        }
      }

      for(i in 1:noAgen){
        if( i == 1){
          Mass[noAgen] <- 0
        }else{
          Mass[i] <- mass[i]*const/sum(mass)
        }
      }


      #calculate force
      if (noAgen==1){
        force <- matrix(0,pop,K)
      }else{
        for(i in 2:iter){
          for(j in 1:K){
            for (l in 2:iter){
              if (U.old[l,j]!=U.old[i,j]){
                force[i,j] <- force[i,j]+ 0.9999 * G * (( Mass[l]* Mass[i] )/ ((U.old[l,j]-U.old[i,j])+0.000001 ))*(U.old[l,j]-U.old[i,j])
              }
            }
          }
        }
      }


      #calculate acceleration
      for(i in 1:iter){
        for (j in 1:K){
          if (Mass[i] != 0){
            A[i,j] <- force[i,j]/Mass[i]
          }
        }
      }


      #calculate velocity
      for(i in 1:iter){
        for (j in 1:K){
          Vm[i,j] <- Velocity[i,j]
          Velocity[i,j] <- Velocity[i,j]+A[i,j]
        }
      }

      #new membership
      for(i in 1:n){
        for(l in 1:K){
          pembilang <- 0
          for(j in 1:p){
            pembilang <- pembilang + abs(data.X[i,j] - V[l,j]) ^ 2
          }
          D[i,l] <- pembilang
        }
      }

      #     distout<-sqrt(dist)
      #     objVal<- membership0*dist


      D <- (D)^(-1/(m-1))
      U.old <- D/rowSums(D)%*%t(matrix(1,K))

      #weigthing membership
      membership <- weighted%*%U.old
      summ <- as.matrix(rowSums(membership))
      #bagi tiap baris matriks U dengan pembagi
      for(j in 1:K){
        for(i in 1:n){
          U.old[i,j] <- (alfa*U.old[i,j])+(beta*(membership[i,j]/summ[i,1]))
        }
      }

      U <- abs(U.old-Velocity)
    }

    if(noAgen > 1){
      if(abs(arrayBestfit[noAgen] - arrayBestfit[noAgen-1]) < threshold){
        break
      }
    }
  }

  #new membership
  for(i in 1:n){
    for(l in 1:K){
      pembilang <- 0
      for(j in 1:p){
        pembilang <- pembilang + abs(data.X[i,j] - V[l,j]) ^ 2
      }
      D[i,l] <- pembilang
    }
  }

  #     distout<-sqrt(dist)
  #     objVal<- membership0*dist


  D <- (D)^(-1/(m-1))
  U <- D/rowSums(D)%*%t(matrix(1,K))
  func.obj <- sum(U ^ m * D)

  func.obj -> func.obj.opt
  U -> U.opt
  V -> V.opt
  D -> D.opt
  ###Labelling###
  colnames(U.opt) = paste("Clust",1:K,sep = " ")
  Clust.desc <- matrix(0,n,p + 1)
  rownames(Clust.desc) <- rownames(X)
  colnames(Clust.desc) <- c(colnames(X),"cluster")
  Clust.desc[,1:p] <- data.X
  for (i in 1:n)
    Clust.desc[i,p + 1] <- which.max(U.opt[i,])
  result <- list()
  result$func.obj <- func.obj.opt
  result$U <- U.opt
  result$V <- V.opt
  result$D <- D.opt
  result$m <- m
  result$data <- data.X
  result$call<-match.call()
  result$Clust.desc <- Clust.desc
  class(result)<-"fgwc-gsa"
  print(result$U)
  result$map <- map
  return(result)
}
