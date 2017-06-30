#' Log Likelihood Function with Constraints Embedded
#'
log.lik.const.1 <- function(pars,
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


  piYW.pop <- c(piYW.pop, 1 - piYW.pop)
  piYM.pop <- c(piYM.pop, 1 - piYM.pop)
  piOW.pop <- c(piOW.pop, 1 - piOW.pop)
  piOM.pop <- c(piOM.pop, 1 - piOM.pop)

  piYW.deg <- c(dnbinom(0:(K - 1), size = piYW.r, prob = 1 - piYW.p),
                1 - sum(dnbinom(0:(K - 1), size = piYW.r, prob = 1 - piYW.p)))
  piYM.deg <- c(dnbinom(0:(K - 1), size = piYM.r, prob = 1 - piYM.p),
                1 - sum(dnbinom(0:(K - 1), size = piYM.r, prob = 1 - piYM.p)))
  piOW.deg <- c(dnbinom(0:(K - 1), size = piOW.r, prob = 1 - piOW.p),
                1 - sum(dnbinom(0:(K - 1), size = piOW.r, prob = 1 - piOW.p)))
  piOM.deg <- c(dnbinom(0:(K - 1), size = piOM.r, prob = 1 - piOM.p),
                1 - sum(dnbinom(0:(K - 1), size = piOM.r, prob = 1 - piOM.p)))

  # Compute intermediate values to make writing the likelihood easier
  muL<-sum(0:K * piYW.deg * piYW.pop[2] * M["YW"],
           0:K * piOW.deg * piOW.pop[2] * M["OW"]) #row totals for les
  muG <- sum(0:K * piYM.deg * piYM.pop[2] * M["YM"],
             0:K * piOM.deg * piOM.pop[2] * M["OM"]) # row totals for gay men
  muBW<-sum(0:K * piYW.deg * piYW.pop[1] * M["YW"],
            0:K * piOW.deg * piOW.pop[1] * M["OW"]) # row totals for bisexual women
  muBM<-sum(0:K * piYM.deg * piYM.pop[1] * M["YM"],
            0:K * piOM.deg * piOM.pop[1] * M["OM"]) # row totals for bisexual men
  mu <-muL + muG + muBW + muBM
  mul<-sum(0:K * piYW.deg * m["HYW"],
           0:K * piOW.deg * m["HOW"]) # scaled row totals for les
  mug<-sum(0:K * piYM.deg * m["HYM"],
           0:K * piOM.deg * m["HOM"]) # scaled row totals for gay
  mubw<-sum(0:K * piYW.deg * m["BYW"],
            0:K * piOW.deg * m["BOW"]) # scaled row totals for bisexual women
  mubm<-sum(0:K * piYM.deg * m["BYM"],
            0:K * piOM.deg * m["BOM"]) # scaled row totals for bisexual men

  # Define Alphas for Convenience
  alphaLG <- alpha[1]
  alphaBWL <- alpha[2]
  alphaBWG <- alpha[3]
  alphaBML <- alpha[4]
  alphaBMG <- alpha[5]
  # Compute additional alphas
  alphaLL <- (muL - alphaLG*muL*muG/mu - alphaBWL*muL*muBW/mu - alphaBML*muL*muBM/mu)*mu/(muL*muL)
  alphaGG <- (muG - alphaLG*muL*muG/mu - alphaBWG*muG*muBW/mu - alphaBMG*muG*muBM/mu)*mu/(muG*muG)
  alphaBWB <- (muBW - alphaBWL*muL*muBW/mu - alphaBWG*muG*muBW/mu)*mu/(muBW*(muBW + muBM))
  alphaBMB <- (muBM - alphaBML*muL*muBM/mu - alphaBMG*muG*muBM/mu)*mu/(muBM*(muBW + muBM))

  # Evaluate the log-likelihood for the given parameter values.  Return the negative log-likelihood
  return(-1 * (
    dmultinom(n[1, 1:3], prob = c(muL / mu * alphaLL,
                                  muG / mu * alphaLG,
                                  (alphaBWL*muBW + alphaBML*muBM) / mu), log = TRUE) +
      dmultinom(n[2, 1:3], prob = c(muL / mu * alphaLG,
                                    muG / mu * alphaGG,
                                    (alphaBWG*muBW + alphaBMG*muBM) / mu), log = TRUE) +
      dmultinom(n[3, 1:3], prob = c(muL / (mu) * alphaBWL,
                                    muG / (mu) * alphaBWG,
                                    (muBW + muBM) * alphaBWB / mu), log = TRUE) +
      dmultinom(n[4, 1:3], prob = c(muL / (mu) * alphaBML,
                                    muG / (mu) * alphaBMG,
                                    (muBW + muBM) * alphaBMB / mu), log = TRUE) +
      dmultinom(DYW, prob = piYW.deg[kYW + 1], log = TRUE) +
      dmultinom(DYM, prob = piYM.deg[kYM + 1], log = TRUE) +
      dmultinom(DOW, prob = piOW.deg[kOW + 1], log = TRUE) +
      dmultinom(DOM, prob = piOM.deg[kOM + 1], log = TRUE) +
      dmultinom(m[grepl("YW", names(m))], prob = piYW.pop, log = TRUE) +
      dmultinom(m[grepl("YM", names(m))], prob = piYM.pop, log = TRUE) +
      dmultinom(m[grepl("OW", names(m))], prob = piOW.pop, log = TRUE) +
      dmultinom(m[grepl("OM", names(m))], prob = piOM.pop, log = TRUE)))
}
