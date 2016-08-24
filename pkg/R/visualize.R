#' Cluster Validity Index
#'
#' @description This function visualize the clustering result
#'
#' @param fgwc result(object) from fgwc clustering
#'
#' @return plot this function will provide visualization using biplot, radar plot and if you use
#' shapefile this function will give maps visualization.
#'
#' @export

visualize <- function(fgwc){


}

function(PC, x="PC1", y="PC2",clust,rowname) {

  # PC being a prcomp object
  data <- data.frame(obsnames=rowname, PC$x)
  clust <- as.character(clust)
  data <- cbind(data,clust)
  # plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      V1 = .7 * mult * (get(x)),
                      V2 = .7 * mult * (get(y))
  )
  ggplot(data, aes(x = PC1, y = PC2),label=TRUE) + geom_point() +
    coord_equal(ratio = 1) +geom_text(data = datapc, aes(x = V1, y = V2, label = varnames), size = 3, vjust = 1) +
    geom_segment(data = datapc, aes(x = 0, y=0, xend = V1, yend = V2),
                 arrow = arrow(length = unit(0.1,"cm")),color = "navy") + geom_text(data = data, aes(group=clust,colour=clust,label = rowname))

}
