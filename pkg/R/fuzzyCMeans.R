#' Fuzzy C-Means
#'
#' @description This function used to perform Fuzzy C-Means of X dataset.
#'
#' @param X data frame n x p
#' @param K specific number of cluster (must be >1)
#' @param m fuzzifier / degree of fuzziness
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
#' @details This function perform Fuzzy C-Means algorithm by Bezdek (1981).
#' Fuzzy C-Means is one of fuzzy clustering methods to clustering dataset
#' become K cluster. Number of cluster (K) must be greater than 1. To control the overlaping
#' or fuzziness of clustering, parameter m must be specified.
#' Maximum iteration and threshold is specific number for convergencing the cluster.
#' Random Number is number that will be used for seeding to firstly generate fuzzy membership matrix.
#' @details Clustering will produce fuzzy membership matrix (U) and fuzzy cluster centroid (V).
#' The greatest value of membership on data point will determine cluster label.
#' Centroid or cluster center can be use to interpret the cluster. Both membership and centroid produced by
#' calculating mathematical distance. Fuzzy C-Means calculate distance with Euclideans norm. So it can be said that cluster
#' will have sperichal shape of geometry.
#'
#' @references Balasko, B., Abonyi, J., & Feil, B. (2002). Fuzzy Clustering and Data Analysis Toolbox: For Use with Matlab. Veszprem, Hungary.
#' @references Gustafson, D. E., & Kessel, W. C. (1978). Fuzzy Clustering With A Fuzzy Covariance Matrix. 761-766.
#' @references Bezdek, J. C., Ehrlich, R., & Full, W. (1984). FCM: The Fuzzy C-Means Clustering Algorithm. Computers and Geosciences Vol 10, 191-203
#'
#' @export
fuzzy.CM<- function(X,K=2,m=2,max.iteration=100,threshold=10^-5,
                     RandomNumber=0) {
  ## Set data
  data.X <- as.matrix(X)
  n <- nrow(data.X)
  p <- ncol(data.X)
  ##Initiation Parameter##
  if (
    (K <= 1) || !(is.numeric(K)) || (K %% ceiling(K) > 0))
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

  U.old <- U + 1
  iteration = 0
  while ((max(abs(U.old - U)) > threshold) &&
         (iteration < max.iteration))
  {
    U.old <- U
    D.old <-D
    V.old<-V
    ## Calculate Centroid
    V <- t(U ^ m) %*% data.X / colSums(U ^ m)
    for (k in 1:K)
    {
      #Distance calculation
      for (i in 1:n)
      {
        D[i,k] = t(data.X[i,] - V[k,]) %*%
          (data.X[i,] -V[k,])
      }
    }
    ##FUZZY PARTITION MATRIX
    for (i in 1:n)
    {
      U[i,] <- 1 /
        (((D[i,]+10^-10) ^ (1 / (m - 1))) *
           sum((1 / (D[i,]+10^-10)) ^ (1 /(m - 1))))
    }
    if(any(is.na(U))==T||any(is.infinite(U))==T)
    {
      U<-U.old
      V<-V.old
      D<-D.old
    }
    for (i in 1:n)
      for (k in 1:K) {
        if (U[i,k] < 0)
          U[i,k] = 0
        else if (U[i,k] > 1)
          U[i,k] = 1
      }
    func.obj = 0
    func.obj = sum(U ^ m * D)
    iteration = iteration + 1
  }
  func.obj -> func.Obj.opt
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
  result$func.obj <- func.Obj.opt
  result$U <- U.opt
  result$V <- V.opt
  result$D <- D.opt
  result$m <- m
  result$call<-match.call()
  result$Clust.desc <- Clust.desc
  class(result)<-"fuzzyclust"
  print(result)
  return(result)
}
