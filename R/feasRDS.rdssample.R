#' Perform Respondent Driven Sampling (RDS)
#'
#' @param net Network to Perform RDS On
#' @param nnodes Number of Nodes in Network to Perform RDS Sample On
#' @param nsamp0 Number of Seeds
#' @param nsamp Goal Sample Size (Is it a problem that we're approaching this as getting a certain size sample?)
#' @param replace Allow Replacement
#' @param coupons Number of Coupons Each Sampled Person Gets
#' @param seed.distribution Probabilities of being selected as seed if not 1/degree (default)
#' @param seed.composition Seed counts by variable
#' @param seed.composition.by Variable differentiating seed counts
#' @param trait.variable "Outcome" of Interest, Needs to be Binary (Diseased = 1)
#' @param recruitment.rates If Not Null, Variable Recruitment Rates for Number of Coupons (Matrix with # of rows = # of "recruitment.rates.by" levels,  2 columns for willingness to spread and willingness to participate named "S" and "C",  row names = levels of "recruitment.rates.by")
#' @param recruitment.rates.by Variable Differentiating Differential Recruitment There is Differential Recruitment
#'
#' @return RDS object
#'
#' @export
#'
feasRDS.rdssample <- function (net, # Network to Perform RDS On
                               nnodes = network.size(net), # Number of Nodes in Network to Perform RDS Sample On
                               nsamp0, # Number of Seeds
                               nsamp, # Goal Sample Size (Is it a problem that we're approaching this as getting a certain size sample?)
                               replace = FALSE, # Allow Replacement
                               coupons, # Number of Coupons Each Sampled Person Gets
                               seed.distribution = NULL, # Probabilities of being selected as seed if not 1/degree (default)
                               seed.composition = NULL, # Seed counts by variable
                               seed.composition.by = NULL, # Variable differentiating seed counts
                               trait.variable = "disease", # "Outcome" of Interest, Needs to be Binary (Diseased = 1)
                               recruitment.rates = NULL, # If Not Null, Variable Recruitment Rates for Number of Coupons
                               # (Matrix with # of rows = # of "recruitment.rates.by" levels,
                               # 2 columns for willingness to spread and willingness to participate named "S" and "C",
                               # row names = levels of "recruitment.rates.by")
                               recruitment.rates.by = NULL) # Variable Differentiating Differential Recruitment There is Differential Recruitment
{
  net <- as.network.uncompressed.rds(net)
  disease <- network::get.vertex.attribute(net, trait.variable)
  gender <- network::get.vertex.attribute(net, "genbin")
  depression <- network::get.vertex.attribute(net, "depbin")
  # Checks That Disease Variable is In Data and Defined
  if (is.null(disease) || all(is.na(disease))) {
    stop(sprintf("No variable called %s appears in the data.",
                 trait.variable))
  }
  # Procedure When Differential Recruitment Rates Are Given
  if (!is.null(recruitment.rates)) {
    if (!is.null(recruitment.rates.by)) {
      rec.var <- network::get.vertex.attribute(net, recruitment.rates.by)
      if (is.null(rec.var) || all(is.na(rec.var))) {
        stop(sprintf("No variable called %s appears in the data.",
                     recruitment.rates.by))
      }
      if (any(sort(unique(rec.var)) != sort(dimnames(recruitment.rates)[[1]]))) {
        stop(sprintf("Levels of %s not the same as row names of recruitment.rates.",
                     recruitment.rates.by))
      }
      if (dim(recruitment.rates)[2] != 2) {
        stop(sprintf("recruitment.rates must have columns equal to 2 for spread and participate"))
      }
    }
  }
  # Get Degrees
  degs <- sapply(net$iel, length) + sapply(net$oel, length)
  nsample <- NULL
  wsample <- NULL
  degsample <- NULL
  dissample <- NULL
  todis <- rep(0, nsamp)
  tonodis <- rep(0, nsamp)
  nominators <- NULL
  # Get Probilities of Seed Selection
  if (is.null(seed.distribution)) {
    # By Degree (Default)
    ps <- sapply(net$iel, length) + sapply(net$oel, length)
  } else {
    ps <- seed.distribution
  }
  if (any(ps < 0)) {
    stop("Negative probabilities of seed selection are not allowed.")
  }
  if (all(ps == 0)) {
    stop("At least one node should have a positive probability of seed selection.\n Are you sure you are selecting some nodes with positive degree?.")
  }

  if (is.null(seed.composition)) {
    if (sum(ps > 0) < nsamp0) {
      nsample <- sample.int(nnodes, size = nsamp0, prob = ps +
                              1e-09, replace = FALSE)
    } else {
      nsample <- sample.int(nnodes, size = nsamp0, prob = ps, replace = FALSE)
    }
  } else {
    seedvar <- network::get.vertex.attribute(net, seed.composition.by)
    nsample <- c()
    vals <- unique(seedvar)[order(unique(seedvar))]
    for (s in 1:length(seed.composition)) {
      if (seed.composition[s] > 0) {
        nsample <- c(nsample,
                     sample(which(seedvar == vals[s]),
                            size = seed.composition[s],
                            prob = ps[seedvar == vals[s]],
                            replace = FALSE))
      }
    }
  }

  wsample <- rep(0, nsamp0)
  nominators <- rep(0, nsamp0)
  refnode <- 1
  while (length(nsample) < nsamp) {
    aaa <- .Call("getNeighborhood_R", net, nsample[refnode],
                 "combined", TRUE, PACKAGE = "network")
    if (!replace) {
      aaa <- setdiff(aaa, nsample)
    }
    prob <- rep(1, length(aaa))

    if (!is.null(recruitment.rates)) {
      if (!is.null(recruitment.rates.by)) {
        # Choose Number of Target Recruits Based on Recruitment Rate of Giver
        target_num_recruits <- sample(c(0, coupons), size = 1,
                                      prob = c(1 - recruitment.rates[rec.var[nsample[refnode]], "S"],
                                               recruitment.rates[rec.var[nsample[refnode]], "S"]))
      } else {
        target_num_recruits <- sample(c(0, coupons), size = 1,
                                      prob = c(1 - recruitment.rates["S"],
                                               recruitment.rates["S"]))
      }
    } else {
      target_num_recruits <- coupons
    }
    if ((nsamp - length(nsample)) >= target_num_recruits) {
      ntosample <- target_num_recruits
    } else {
      ntosample <- nsamp - length(nsample)
    }
    if (ntosample > length(aaa)) {
      ntosample <- length(aaa)
    }
    if (length(aaa) > 1) {
      thissamp <- sample(aaa, size = ntosample, prob = prob)
      if (!is.null(recruitment.rates)) {
        if (!is.null(recruitment.rates.by)) {
          thissamp <- thissamp[runif(ntosample, 0, 1) < recruitment.rates[rec.var[nsample[refnode]], "C"]]
        } else {
          thissamp <- thissamp[runif(ntosample, 0, 1) < recruitment.rates["C"]]
        }
      }
      ntosample <- length(thissamp)
      if (is.na(thissamp[1])) {
        thissamp <- NULL
        ntosample <- 0
      }
    } else {
      thissamp <- NULL
    }
    if (length(aaa) == 1) {
      thissamp <- aaa
      if (!is.null(recruitment.rates)) {
        if (!is.null(recruitment.rates.by)) {
          thissamp <- thissamp[runif(1, 0, 1) < recruitment.rates[rec.var[nsample[refnode]], "C"]]
        } else {
          thissamp <- thissamp[runif(1, 0, 1) < recruitment.rates["C"]]
        }
      }
      ntosample <- length(ntosample)
      if (is.na(thissamp[1])) {
        thissamp <- NULL
        ntosample <- 0
      }
    }
    if (ntosample > 0) {
      nsample <- c(nsample, thissamp)
      wsample <- c(wsample, rep((wsample[refnode] + 1),
                                ntosample))
      nominators <- c(nominators, rep(nsample[refnode],
                                      ntosample))
      for (nn in 1:ntosample) {
        if (disease[thissamp[nn]] == 1) {
          todis[refnode] <- todis[refnode] + 1
        }
        if (disease[thissamp[nn]] == 0) {
          tonodis[refnode] <- tonodis[refnode] + 1
        }
      }
    }
    refnode <- refnode + 1
    # This Means Sampling Has Gotten Stuck
    if (refnode > length(nsample)) {
      cat(sprintf("trapped, sampled so far=%d.\n",
                  length(nsample)))
      break

    }
  }
  # Collect Information About Sample
  degsample <- degs[nsample]
  dissample <- disease[nsample]
  gensample <- gender[nsample]
  depsample <- depression[nsample]
  if (!is.null(recruitment.rates.by)) {
    recsample <- rec.var[nsample]
  } else {
    recsample <- rep(NA, length(degsample))
  }
  list(nsample = nsample, wsample = wsample, nominators = nominators,
         degsample = degsample, dissample = dissample, todis = todis,
         tonodis = tonodis, recsample = recsample, gensample = gensample,
         depsample = depsample)
}
