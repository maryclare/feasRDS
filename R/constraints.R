#' Defines a function that evaluates the non-constant part of a constraint
#' at the provided parameter and data values
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
#' @return A vector of values corresponding to the nonconstant parts of the constraints evaluated at the provided data and parameter values
#'
#' @export
constraints <- function(pars,
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

  # Constraints that ensure that alpha parameters ensure consistency between expected mixing totals and marginals of the mixing matrix
  z1 <- muL -
    muL * muL / mu * alpha[1] -
    muL * muG / mu * alpha[2] -
    muL * muBW / mu * alpha[4] -
    muL * muBM / mu * alpha[7]
  z2 <- muG -
    muG * muL / mu * alpha[2] -
    muG * muG / mu * alpha[3] -
    muG * muBW / mu * alpha[5] -
    muG * muBM / mu * alpha[8]
  z3 <- muBW -
    muL * (alpha[4]*muBW) / mu -
    muG * (alpha[5]*muBW) / mu -
    muBW * (muBW + muBM) / mu * alpha[6]
  z4 <- muBM -
    muL * (alpha[7]*muBM) / mu -
    muG * (alpha[8]*muBM) / mu -
    muBM * (muBW + muBM) / mu * alpha[9]

  # Ensure multinomial probabilities sum to 1
  z5 <- sum(piYW.pop)
  z6 <- sum(piYM.pop)
  z7 <- sum(piOW.pop)
  z8 <- sum(piOM.pop)

  # Return the constraint values
  return(c(z1, z2, z3, z4, z5, z6, z7, z8))
}
