#' Performs Reed-Molloy algorithm to initialize degrees
#'
#' @param degrees Vector of desired degree counts to be matched in the simulated network, where each element is the degree count of a synthetic subject
#'
#' @return A network object with a degree distribution equal or close to the desired degree distribution
#'
#' @export
reedm <- function(degrees){

  od <- (degrees>0)
  degs <- degrees[od]

  done<-"no"
  while(done!="yes"){  #want to use reedmolloy to find a suitable network.  iterate finding of counts at the same time.

    cat(sprintf("Trial: Is stubs even? The total is %d\n",sum(degs)))
    #
    if(sum(degs)%%2|done=="bad"){
      #   degs <- roundstoc.RDS(degs+stats::runif(length(degs),-0.2,0.2))
      d <- sample(which(degs > 1),size=1)
      degs[d] <- degs[d] + sample(c(-1,1),size=1)
      done<-"no"
    }else{
      done <- "yes"
      sm0 <- try(igraph::get.edgelist(igraph::degree.sequence.game(degs,method="vl")))
      if(inherits(sm0,"try-error")){ done<-"bad" }else{
        #      Add back 0 degrees
        a=cumsum(!od)[od]+seq_along(degs)
        net <- network.initialize(length(degrees), directed=FALSE)
        net <- network::add.edges(x=net,
                                  tail=as.list(a[sm0[,1]]),
                                  head=as.list(a[sm0[,2]]))
      }
    }
  }
  net
}
