#' Data Scalling
#'
#' @description Provide data scalling using z-transform, zero to one scalling and minus one to one scalling
#'
#' @param data matrix data
#' @param method scalling technique use "z" for z-transform, "zerotoone" for zero to one scalling and "oneminuseone" minus one to one scalling
#'
#' @return scalled matrix data
#'
#' @seealso \code{\link{fgwc}} for standard Fuzzy Geographically Weighted Clustering,
#' \code{\link{fgwc.gsa}} for optimize using Gravitational Search Algorithm,
#' \code{\link{spClustIndex}} for cluser validation,
#' \code{\link{visualize}} for cluster visualizatiion
#'
#' @examples
#' #load data
#' data <- example
#'
#' #zero to one scalling
#' data <- scale(data)
#' data <- scale(data,method="zerotoone")
#'
#' #z-transform
#' data <- scale(data,method="z")
#'
#' #minus one to one scalling
#' data <- scale(data,method="oneminuseone")
#'
#' @export
#'

scale <- function(data,method="zerotoone"){

  if(method=="z"){
    data <- scale.default(data)
    return (data)
  }

  if(method=="zerotoone"){
    data <- normalizeZeroOne(data)
    return (data)
  }

  if(method=="oneminuseone"){
    data <- normalizeOneMinOne(data)
    return (data)
  }

}

normalizeOneMinOne <- function(data){

  for(j in 1:ncol(data)){
    x <- data[,j]
    mx <- max(x)
    mn <- min(x)

    scl1 <- (mx-mn)/2
    scl2 <- (mx+mn)/2
    for(i in 1:length(x)){
      data[i,j] <- (data[i,j]-scl2)/scl1
    }
  }


  return(data)
}

normalizeZeroOne <- function(data){
  for(j in 1:ncol(data)){
    x <- data[,j]
    mx <- max(x)
    mn <- min(x)

    scl1 <- (mx-mn)
    for(i in 1:length(x)){
      data[i,j] <- (data[i,j]-mn)/scl1
    }
  }

  return(data)
}
