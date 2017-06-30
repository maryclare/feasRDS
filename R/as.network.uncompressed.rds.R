#' Function that extracts network features
#'
as.network.uncompressed.rds<-function(x,
                                      na.rm=FALSE, edge.check=FALSE){
  #Initialize the network object
  if(class(x)=="network"){return(x)}
  if(is.null(attr(x,"vnames"))){
    warning("as.network.uncompressed.rds input must be a compressed network, or a network.\n Returning the original object.\n")
    return(x)
  }
  n<-attr(x,"n")
  directed<-attr(x,"directed")
  g<-network.initialize(n,directed=directed)
  #Call the specific coercion routine, depending on matrix type
  # g<-network.edgelist(x,g,na.rm=na.rm,edge.check=edge.check)
  g<-network::add.edges(g,as.list(x[,1]),as.list(x[,2]),edge.check=edge.check)
  va <- attr(x,"vertex.attributes")
  if(length(va)>0){
    for (i in (1:length(va))){
      g <- network::set.vertex.attribute(g,names(va)[i], va[[i]])
    }
  }
  #Return the result
  g
}
