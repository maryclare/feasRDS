#' Function that takes data for one MSA and returns degree distributions, mixing matrices, starting values...
#'
#' @param data.sub Data for one MSA
#' @param K Degree cutoff to be used
#'
#' @return List of objects derived from the data used in maximum likelihood estimation and/or RDS simulation
#'
#' @export
get.data <- function(data.sub = data.sub, # Appropriate subset of the data
                     K = 45 # Degree cutoff to be used
) {

  ### Get sample sizes
  # Initialize matrix of sample sizes by
  # sexual identity, gender and age group
  m.0 <- expand.grid("bisexual" = unique(data.sub$bisexual),
                     "gender4" = unique(data.sub$gender4),
                     "ageimp.ind" = unique(data.sub$ageimp.ind))
  # Get sample sizes from data
  m <- aggregate(data.sub$id, list(bisexual = data.sub$bisexual,
                                   gender4 = data.sub$gender4,
                                   ageimp.ind = data.sub$ageimp.ind),
                 function(x) {length(x)})
  names(m)[ncol(m)] <- "m"
  # Merge sample sizes from data into initalized matrix of
  # sample sizes
  m <- merge(m, m.0, all = TRUE)
  # Fill in any unrepresented groups with 0's
  m$m[is.na(m$m)] <- 0
  # Put sample sizes in desired order
  m <- m[order(m$ageimp.ind, m$gender4, m$bisexual), ]
  m <- m$m
  # Name the vector of sample sizes
  names(m) <- c("BYW", "HYW",
                "BYM", "HYM",
                "BOW", "HOW",
                "BOM", "HOM")

  M <- c(0.25, 0.45, 0.15, 0.15)*5000

  names(M) <- c("YW",
                "YM",
                "OW",
                "OM")

  # Get sample degree distributions, mixing totals

  # First, need to truncate degree counts, rescale
  data.sub$les50.imp.round.trunc <- data.sub$les50.imp.round
  data.sub$gay50.imp.round.trunc <- data.sub$gay50.imp.round
  data.sub$bi50.imp.round.trunc <- data.sub$bi50.imp.round

  # Calculate degree from stratified alter counts
  data.sub$tot50.imp.round <- data.sub$les50.imp.round + data.sub$gay50.imp.round + data.sub$bi50.imp.round

  # Rescale
  for (v in c("les", "gay", "bi")) {
    data.sub[data.sub$tot50.imp.round > K, paste(v, "50.imp.round.trunc", sep = "")] <- (K/data.sub[data.sub$tot50.imp.round > K, "tot50.imp.round"])*data.sub[data.sub$tot50.imp.round > K, paste(v, "50.imp.round", sep = "")]
  }
  # Truncate
  data.sub$tot50.imp.round.trunc <- data.sub$tot50.imp.round
  data.sub$tot50.imp.round.trunc[data.sub$tot50.imp.round.trunc > K] <- K

  # Get vectors of degree values that are represented
  kYW <- as.numeric(names(table(data.sub$tot50.imp.round.trunc[data.sub$gender4 == "female" & data.sub$ageimp.ind == 0])))
  kYM <- as.numeric(names(table(data.sub$tot50.imp.round.trunc[data.sub$gender4 == "male" & data.sub$ageimp.ind == 0])))
  kOW <- as.numeric(names(table(data.sub$tot50.imp.round.trunc[data.sub$gender4 == "female" & data.sub$ageimp.ind == 1])))
  kOM <- as.numeric(names(table(data.sub$tot50.imp.round.trunc[data.sub$gender4 == "male" & data.sub$ageimp.ind == 1])))

  # Get vectors of counts for each degree value
  DYW <- as.vector(table(data.sub$tot50.imp.round.trunc[data.sub$gender4 == "female" & data.sub$ageimp.ind == 0]))
  DYM <- as.vector(table(data.sub$tot50.imp.round.trunc[data.sub$gender4 == "male" & data.sub$ageimp.ind == 0]))
  DOW <- as.vector(table(data.sub$tot50.imp.round.trunc[data.sub$gender4 == "female" & data.sub$ageimp.ind == 1]))
  DOM <- as.vector(table(data.sub$tot50.imp.round.trunc[data.sub$gender4 == "male" & data.sub$ageimp.ind == 1]))

  # Set names of counts to degree values
  names(DYW) <- kYW
  names(DYM) <- kYM
  names(DOW) <- kOW
  names(DOM) <- kOM

  # Put sample degree distributions in a list
  D <- list(DYW, DYM, DOW, DOM)
  names(D) <- c("YW", "YM", "OW", "OM")

  # Get matrix of sample mixing totals
  n <- round(cbind(aggregate(data.sub$les50.imp.round.trunc, list(gender4 = data.sub$gender4,
                                                                  bisexual = data.sub$bisexual), sum)$x,
                   aggregate(data.sub$gay50.imp.round.trunc, list(gender4 = data.sub$gender4,
                                                                  bisexual = data.sub$bisexual), sum)$x,
                   aggregate(data.sub$bi50.imp.round.trunc, list(gender4 = data.sub$gender4,
                                                                 bisexual = data.sub$bisexual), sum)$x)[c(3, 4, 1, 2), ])

  row.names(n) <- c("L", "G", "BW", "BM")
  colnames(n) <- c("L", "G", "B")

  ### Get starting values
  # Take vectors of degree reports that have not been imputed
  # to set negative binomial parameters
  dYW <- data.sub$tot50.imp.round[data.sub$gender4 == "female" & data.sub$ageimp.ind == 0]
  dYM <- data.sub$tot50.imp.round[data.sub$gender4 == "male" & data.sub$ageimp.ind == 0]
  dOW <- data.sub$tot50.imp.round[data.sub$gender4 == "female" & data.sub$ageimp.ind == 1]
  dOM <- data.sub$tot50.imp.round[data.sub$gender4 == "male" & data.sub$ageimp.ind == 1]

  # Get Moment Estimators for Negative Binomial Distribution
  mYW <- mean(dYW)
  vYW <- var(dYW)
  pYW <- -1*(mYW/vYW - 1)
  rYW <- mYW*(1 - pYW)/pYW
  mOW <- mean(dOW)
  vOW <- var(dOW)
  pOW <- -1*(mOW/vOW - 1)
  rOW <- mOW*(1 - pOW)/pOW
  mYM <- mean(dYM)
  vYM <- var(dYM)
  pYM <- -1*(mYM/vYM - 1)
  rYM <- mYM*(1 - pYM)/pYM
  mOM <- mean(dOM)
  vOM <- var(dOM)
  pOM <- -1*(mOM/vOM - 1)
  rOM <- mOM*(1 - pOM)/pOM

  # Get Moment Estimators for Parameters
  start.values.nb  <- c(rYW, pYW,
                        rYM, pYM,
                        rOW, pOW,
                        rOM, pOM,
                        rep(1, 9),
                        (m + 1)/c(rep(sum((m[grepl("YW", names(m))] + 1)), 2),
                                  rep(sum((m[grepl("YM", names(m))] + 1)), 2),
                                  rep(sum((m[grepl("OW", names(m))] + 1)), 2),
                                  rep(sum((m[grepl("OM", names(m))] + 1)), 2)))

  names(start.values.nb) <- c(c(paste("piYW", c("r", "p"), sep = "."),
                                paste("piYM", c("r", "p"), sep = "."),
                                paste("piOW", c("r", "p"), sep = "."),
                                paste("piOM", c("r", "p"), sep = ".")),
                              paste("alpha", c("LL", "LG", "GG", "BWL", "BWG", "BWB", "BML", "BMG", "BMB"), sep = ""),
                              paste("pi", c("BYW", "HYW",
                                            "BYM", "HYM",
                                            "BOW", "HOW",
                                            "BOM", "HOM"), sep = ""))

  ### Get recruitment rates for RDS simulation by gender,
  ### sexual identity and age group
  # Make variables into indicators
  data.sub$spread.ind <- data.sub$spread == "willing to spread word about future project"
  data.sub$contact.ind <- data.sub$contact == "willing to be contacted in future"

  # Get aggregate rates of willingness to participate and recruit
  spread <- aggregate(data.sub$spread.ind, list(gender4 = data.sub$gender4,
                                                bisexual = data.sub$bisexual,
                                                ageimp.ind = data.sub$ageimp.ind),
                      mean, na.rm = TRUE)
  contact <- aggregate(data.sub$contact.ind, list(gender4 = data.sub$gender4,
                                                  bisexual = data.sub$bisexual,
                                                  ageimp.ind = data.sub$ageimp.ind),
                       mean, na.rm = TRUE)

  # Combine rates into one data frame and name rows and columns
  recruit <- merge(spread, contact, by = c("gender4", "bisexual", "ageimp.ind"))
  row.names(recruit)[recruit$gender4 == "female" &
                       recruit$bisexual == "bisexual" &
                       recruit$ageimp.ind == FALSE] <- "BYW"
  row.names(recruit)[recruit$gender4 == "male" &
                       recruit$bisexual == "bisexual" &
                       recruit$ageimp.ind == FALSE] <- "BYM"
  row.names(recruit)[recruit$gender4 == "female" &
                       recruit$bisexual == "not bisexual" &
                       recruit$ageimp.ind == FALSE] <- "HYW"
  row.names(recruit)[recruit$gender4 == "male" &
                       recruit$bisexual == "not bisexual" &
                       recruit$ageimp.ind == FALSE] <- "HYM"
  row.names(recruit)[recruit$gender4 == "female" &
                       recruit$bisexual == "bisexual" &
                       recruit$ageimp.ind == TRUE] <- "BOW"
  row.names(recruit)[recruit$gender4 == "male" &
                       recruit$bisexual == "bisexual" &
                       recruit$ageimp.ind == TRUE] <- "BOM"
  row.names(recruit)[recruit$gender4 == "female" &
                       recruit$bisexual == "not bisexual" &
                       recruit$ageimp.ind == TRUE] <- "HOW"
  row.names(recruit)[recruit$gender4 == "male" &
                       recruit$bisexual == "not bisexual" &
                       recruit$ageimp.ind == TRUE] <- "HOM"
  recruit <- recruit[, c("x.x", "x.y")]
  colnames(recruit) <- c("S", "C")
  recruit.pool <- c(mean(data.sub$spread.ind, na.rm = TRUE), mean(data.sub$contact.ind, na.rm = TRUE))
  names(recruit.pool) <- c("S", "C")

  # Get depression rates by age, gender and degree for
  # assignment of depression status
  # Simple GLM of depression
  depreg <- glm(depressed ~ ageimp.ind*factor(gender4)*I(log(tot50.imp.round + 1)), data = data.sub,
                family=binomial(link=logit))
  # Get predicted depression rates
  deppred <- predict(depreg,
                     data.frame(ageimp.ind = rep(c(0, 1, 0, 1), each = length(0:K)),
                                gender4 = rep(c("female", "female", "male", "male"), each = length(0:K)),
                                tot50.imp.round = rep(0:K, times = 4)))
  deprate <- exp(deppred)/(1 + exp(deppred))
  # Get depression rates
  names(deprate) <- paste(rep(0:K, times = 4),
                          rep(c("YW", "OW", "YM", "OM"), each = length(0:K)),
                          sep = ".")

  return(list("m" = m,
              "M" = M,
              "D" = D,
              "n" = n,
              "start" = start.values.nb,
              "part.ref" = recruit,
              "dep.rate" = deprate))
}

