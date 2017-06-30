#' Function to get standard errors for parameters estimated by maximizing "loglik" function
#'
#' @param pars Vector of estimates of the unknown parameters, vector must have names(pars) that identify the individual elements consistent with  names(work[["start"]]), where the object work is output of the get.data function
#'
#' @return Vector of standard error estimates for each parameter in the "pars" vector
#'
#' @export
get.ses <- function(pars) {

  # Need to break up parameters into three groups to embed constraints in likelihood
  pars.se.1 <- pars[c((1:length(pars))[grepl("pi", names(pars)) & !grepl("H", names(pars)) & !grepl("B", names(pars))],
                      which(names(pars) == "alphaLG"),
                      which(names(pars) == "alphaBWL"),
                      which(names(pars) == "alphaBWG"),
                      which(names(pars) == "alphaBML"),
                      which(names(pars) == "alphaBMG"),
                      (1:length(pars))[grepl("piB", names(pars))])]
  pars.se.2 <- pars[c((1:length(pars))[grepl("pi", names(pars)) & !grepl("H", names(pars)) & !grepl("B", names(pars))],
                      which(names(pars) == "alphaLL"),
                      which(names(pars) == "alphaGG"),
                      which(names(pars) == "alphaBWL"),
                      which(names(pars) == "alphaBML"),
                      which(names(pars) == "alphaBMB"),
                      (1:length(pars))[grepl("piH", names(pars))])]
  pars.se.3 <- pars[c((1:length(pars))[grepl("pi", names(pars)) & !grepl("H", names(pars)) & !grepl("B", names(pars))],
                      which(names(pars) == "alphaLL"),
                      which(names(pars) == "alphaGG"),
                      which(names(pars) == "alphaBWG"),
                      which(names(pars) == "alphaBWB"),
                      which(names(pars) == "alphaBMG"),
                      (1:length(pars))[grepl("piB", names(pars))])]

  # Numerically compute hessian for each group of parameters
  hess.1 <- optimHess(pars.se.1,
                      fn = log.lik.const.1,
                      M = work[["M"]],
                      m = work[["m"]],
                      n = work[["n"]],
                      kYW = as.numeric(names(work[["D"]][["YW"]])),
                      kYM = as.numeric(names(work[["D"]][["YM"]])),
                      kOW = as.numeric(names(work[["D"]][["OW"]])),
                      kOM = as.numeric(names(work[["D"]][["OM"]])),
                      DYW = work[["D"]][["YW"]],
                      DYM = work[["D"]][["YM"]],
                      DOW = work[["D"]][["OW"]],
                      DOM = work[["D"]][["OM"]],
                      K = K,
                      control=list(ndeps= pars.se.1 / 10000))

  hess.2 <- optimHess(pars.se.2,
                      fn = log.lik.const.2,
                      M = work[["M"]],
                      m = work[["m"]],
                      n = work[["n"]],
                      kYW = as.numeric(names(work[["D"]][["YW"]])),
                      kYM = as.numeric(names(work[["D"]][["YM"]])),
                      kOW = as.numeric(names(work[["D"]][["OW"]])),
                      kOM = as.numeric(names(work[["D"]][["OM"]])),
                      DYW = work[["D"]][["YW"]],
                      DYM = work[["D"]][["YM"]],
                      DOW = work[["D"]][["OW"]],
                      DOM = work[["D"]][["OM"]],
                      K = K,
                      control=list(ndeps= pars.se.2 / 10000))

  hess.3 <- optimHess(pars.se.3, fn = log.lik.const.3,
                      M = work[["M"]],
                      m = work[["m"]],
                      n = work[["n"]],
                      kYW = as.numeric(names(work[["D"]][["YW"]])),
                      kYM = as.numeric(names(work[["D"]][["YM"]])),
                      kOW = as.numeric(names(work[["D"]][["OW"]])),
                      kOM = as.numeric(names(work[["D"]][["OM"]])),
                      DYW = work[["D"]][["YW"]],
                      DYM = work[["D"]][["YM"]],
                      DOW = work[["D"]][["OW"]],
                      DOM = work[["D"]][["OM"]],
                      K = K,
                      control=list(ndeps= pars.se.3 / 10000))

  # Things that should be the same look the same
  se.1 <- sqrt(diag(solve(hess.1)))
  se.2 <- sqrt(diag(solve(hess.2)))
  se.3 <- sqrt(diag(solve(hess.3)))

  # Paste all the parameters together
  se <- c(se.1[grepl("piYW", names(se.1))],
          se.2[grepl("piYW", names(se.2)) & !names(se.2) %in% names(se.1)],
          se.1[grepl("piYM", names(se.1))],
          se.2[grepl("piYM", names(se.2)) & !names(se.2) %in% names(se.1)],
          se.1[grepl("piOW", names(se.1))],
          se.2[grepl("piOW", names(se.2)) & !names(se.2) %in% names(se.1)],
          se.1[grepl("piOM", names(se.1))],
          se.2[grepl("piOM", names(se.2)) & !names(se.2) %in% names(se.1)])

  pi.se <- se[grepl("pi", names(se)) & !grepl("H", names(se)) & !grepl("B", names(se))]
  pop.se <- c(se.1[grepl("piBYW", names(se.1))],
              se.2[grepl("piHYW", names(se.2))],
              se.1[grepl("piBYM", names(se.1))],
              se.2[grepl("piHYM", names(se.2))],
              se.1[grepl("piBOW", names(se.1))],
              se.2[grepl("piHOW", names(se.2))],
              se.1[grepl("piBOM", names(se.1))],
              se.2[grepl("piHOM", names(se.2))])
  alph.se <- c(c(se.1, se.2, se.3)[grepl("alphaLL", names(c(se.1, se.2, se.3)))][1],
               c(se.1, se.2, se.3)[grepl("alphaLG", names(c(se.1, se.2, se.3)))][1],
               c(se.1, se.2, se.3)[grepl("alphaGG", names(c(se.1, se.2, se.3)))][1],
               c(se.1, se.2, se.3)[grepl("alphaBWL", names(c(se.1, se.2, se.3)))][1],
               c(se.1, se.2, se.3)[grepl("alphaBWG", names(c(se.1, se.2, se.3)))][1],
               c(se.1, se.2, se.3)[grepl("alphaBWB", names(c(se.1, se.2, se.3)))][1],
               c(se.1, se.2, se.3)[grepl("alphaBML", names(c(se.1, se.2, se.3)))][1],
               c(se.1, se.2, se.3)[grepl("alphaBMG", names(c(se.1, se.2, se.3)))][1],
               c(se.1, se.2, se.3)[grepl("alphaBMB", names(c(se.1, se.2, se.3)))][1])

  return(list("se" = se,
              "alph.se" = alph.se,
              "pi.se" = pi.se,
              "pop.se" = pop.se))

}
