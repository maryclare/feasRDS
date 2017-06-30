#' Function that takes object produced by feasRDS.rdssample and produces a nice data frame of RDS data
#'
#' @param x Object returned by rdssample.mc
#' @param population.size Optional value equal to population size if supplied
#' @return A data frame of RDS data
#' @export
feasRDS.to.data.frame<-function(x, population.size=NULL){
  as.rds.data.frame(
    df=data.frame(wave=x$wsample,degree=x$degsample,bibin=x$dissample, type = x$recsample,
                  id=x$nsample,recruiter.id=x$nominators,
                  genbin = x$gensample, depbin = x$depsample),
    population.size=population.size,network.size="degree"
  )
}
