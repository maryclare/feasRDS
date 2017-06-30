#' Returns the log likelihood at the provided parameter and data values
#'
#' @param pars Vector of starting values for estimating unknown parameters, vector must have names(pars) that identify the individual elements consistent with  names(work[["start"]]), where the object work is output of the get.data function
#' @param m Vector of sample sizes by gender, sexual identity and age group, vector must have names(m) that identify the individual elements consistent with  names(work[["m"]]), where the object work is output of the get.data function
#' @param M Vector of sample sizes by gender and age group, vector must have names(M) that identify the individual elements consistent with  names(work[["M"]]), where the object work is output of the get.data function
#' @param n 4 by 3 matrix of sample mixing totals, matrix must have row.names(n) and column.names(n) that identify the individual elements consistent with  names(work[["n"]]), where the object work is output of the get.data function
#' @param kYW Vector of observed degree values in the sample for younger women, ordered from lowest to highest
#' @param kYM Vector of observed degree values in the sample for younger men, ordered from lowest to highest
#' @param kOW Vector of observed degree values in the sample for older women, ordered from lowest to highest
#' @param kOM Vector of observed degree values in the sample for older men, ordered from lowest to highest
#' @param DYW Vector of observed degree counts in the sample for younger women, ordered in increasing order of corresponding degree value
#' @param DYM Vector of observed degree counts in the sample for younger men, ordered in increasing order of corresponding degree value
#' @param DOW Vector of observed degree counts in the sample for older women, ordered in increasing order of corresponding degree value
#' @param DOM Vector of observed degree counts in the sample for older men, ordered in increasing order of corresponding degree value
#' @param K Degree value cutoff
#'
#' @return The log likelihood evaluated at the provided data and parameter values
#'
#' @export
loglik <- function(pars,
                           m,
                           M,
                           n,
                           kYW, kYM, kOW, kOM,
                           DYW, DYM, DOW, DOM,
                           K)
{

  # Extract degree distribution parameters and cell-specific dependence parameters
  if (!is.null(names(pars))) {
    piYW.r <- pars["piYW.r"]
    piYM.r <- pars["piYM.r"]
    piOW.r <- pars["piOW.r"]
    piOM.r <- pars["piOM.r"]
    piYW.p <- pars["piYW.p"]
    piYM.p <- pars["piYM.p"]
    piOW.p <- pars["piOW.p"]
    piOM.p <- pars["piOM.p"]
    alpha <- pars[grepl("alpha", names(pars))]
    piYW.pop <- pars[grepl("BYW", names(pars)) | grepl("HYW", names(pars))]
    piYM.pop <- pars[grepl("BYM", names(pars)) | grepl("HYM", names(pars))]
    piOW.pop <- pars[grepl("BOW", names(pars)) | grepl("HOW", names(pars))]
    piOM.pop <- pars[grepl("BOM", names(pars)) | grepl("HOM", names(pars))]
  }

  piYW.deg <- c(dnbinom(0:(K - 1), size = piYW.r, prob = 1 - piYW.p),
                1 - sum(dnbinom(0:(K - 1), size = piYW.r, prob = 1 - piYW.p)))
  piYM.deg <- c(dnbinom(0:(K - 1), size = piYM.r, prob = 1 - piYM.p),
                1 - sum(dnbinom(0:(K - 1), size = piYM.r, prob = 1 - piYM.p)))
  piOW.deg <- c(dnbinom(0:(K - 1), size = piOW.r, prob = 1 - piOW.p),
                1 - sum(dnbinom(0:(K - 1), size = piOW.r, prob = 1 - piOW.p)))
  piOM.deg <- c(dnbinom(0:(K - 1), size = piOM.r, prob = 1 - piOM.p),
                1 - sum(dnbinom(0:(K - 1), size = piOM.r, prob = 1 - piOM.p)))

  # Compute intermediate values to make writing the likelihood easier
  muL<-sum(0:K * piYW.deg * piYW.pop["piHYW"] * M["YW"],
           0:K * piOW.deg * piOW.pop["piHOW"] * M["OW"]) #row totals for les
  muG <- sum(0:K * piYM.deg * piYM.pop["piHYM"] * M["YM"],
             0:K * piOM.deg * piOM.pop["piHOM"] * M["OM"]) # row totals for gay men
  muBW<-sum(0:K * piYW.deg * piYW.pop["piBYW"] * M["YW"],
            0:K * piOW.deg * piOW.pop["piBOW"] * M["OW"]) # row totals for bisexual women
  muBM<-sum(0:K * piYM.deg * piYM.pop["piBYM"] * M["YM"],
            0:K * piOM.deg * piOM.pop["piBOM"] * M["OM"]) # row totals for bisexual men
  mu <-muL + muG + muBW + muBM
  mul<-sum(0:K * piYW.deg * m["HYW"],
           0:K * piOW.deg * m["HOW"]) # scaled row totals for les
  mug<-sum(0:K * piYM.deg * m["HYM"],
           0:K * piOM.deg * m["HOM"]) # scaled row totals for gay
  mubw<-sum(0:K * piYW.deg * m["BYW"],
            0:K * piOW.deg * m["BOW"]) # scaled row totals for bisexual women
  mubm<-sum(0:K * piYM.deg * m["BYM"],
            0:K * piOM.deg * m["BOM"]) # scaled row totals for bisexual men

  # Evaluate the log-likelihood for the given parameter values.  Return the negative log-likelihood
  return(-1 * (
    dmultinom(n[1, 1:3], prob = c(muL / mu * alpha[1],
                                  muG / mu * alpha[2],
                                  (alpha[4]*muBW + alpha[7]*muBM) / mu), log = TRUE) +
      dmultinom(n[2, 1:3], prob = c(muL / mu * alpha[2],
                                    muG / mu * alpha[3],
                                    (alpha[5]*muBW + alpha[8]*muBM) / mu), log = TRUE) +
      dmultinom(n[3, 1:3], prob = c(muL / (mu) * alpha[4],
                                    muG / (mu) * alpha[5],
                                    (muBW + muBM) * alpha[6] / mu), log = TRUE) +
      dmultinom(n[4, 1:3], prob = c(muL / (mu) * alpha[7],
                                    muG / (mu) * alpha[8],
                                    (muBW + muBM) * alpha[9] / mu), log = TRUE) +
      dmultinom(DYW, prob = piYW.deg[kYW + 1], log = TRUE) +
      dmultinom(DYM, prob = piYM.deg[kYM + 1], log = TRUE) +
      dmultinom(DOW, prob = piOW.deg[kOW + 1], log = TRUE) +
      dmultinom(DOM, prob = piOM.deg[kOM + 1], log = TRUE) +
      dmultinom(m[grepl("YW", names(m))], prob = piYW.pop, log = TRUE) +
      dmultinom(m[grepl("YM", names(m))], prob = piYM.pop, log = TRUE) +
      dmultinom(m[grepl("OW", names(m))], prob = piOW.pop, log = TRUE) +
      dmultinom(m[grepl("OM", names(m))], prob = piOM.pop, log = TRUE)))
}
