# this script contains functions that compute posterior-predictive summary statistics
# under a constant or interval-specific geographical model
# specifically, two such statistics are focused here, including the time-slice tipwise multinomial statistic and time-slice parsimony statistic
# it also contains functions that can reconstruct parsimony history (which is in turn used to compute time-slice parsimony statistic), 
# as well as a function that computes posterior-predictive p-value based on the posterior-predictive distribution of a summary statistic
library(phangorn)
library(stringr)
library(pbapply)

#' Compute the piecewise constant tipwise multinomial statistics for one interval
#' @param nvec a vector of states of all tips
#' @param mvec a vector of states of the tips in a given time interval
#' @return the piecewise constant tipwise multinomial statistic of this interval
multinomiallikelihood_tipwise_calculator <- function(nvec, mvec = NULL) {
  
  nvec <- nvec[nvec != "?"]
  if (is.null(mvec)) { # no interval-specific tips, then all the tips fall in one interval
    return(sum(table(nvec) * log(table(nvec))) - length(nvec) * log(length(nvec)))
  } else {
    if (length(mvec) > 0) {
      mvec <- mvec[mvec != "?"]
      pvec <- table(nvec) / length(nvec)
      nmvec <- table(mvec)
      return(sum(nmvec * log(pvec[names(nmvec)])))
    } else {
      return(0)
    }
  }
}

#' Compute the piecewise constant tipwise multinomial statistics
#' @param tree a bifuracting tree of class phylo
#' @param tipstates the observed tip states
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model
#' @return a vector (each element correspond to a time interval, from present to ancient) of piecewise constant tipwise multinomial statistics
multinomiallikelihood_tipwise_epoch_calculator <- function(tree, tipstates, Q_ages) {
  
  tips <- tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]]
  node_times <- ape::node.depth.edgelength(tree)
  node_ages <- max(node_times) - node_times
  
  Q_times <- max(node_ages) - sort(Q_ages, decreasing = T)
  epoch_num <- length(Q_times) + 1L
  tipstates <- tipstates[tree$tip.label]
  
  tipepoch_idx <- findInterval(node_times[tips], Q_times) + 1L
  # fetch the tips in each interval and perform the calculation
  return(rev(sapply(1:epoch_num, function(m) multinomiallikelihood_tipwise_calculator(nvec = tipstates, mvec = tipstates[tipepoch_idx == m]))))
}

#' Reconstruct a parsimony history of a discrete-character over a bifurcation tree conditioning on the observed tip states and the simulated tip states, respectively
#' @param trees a list of bifuracting trees (each of them is a phylo object with an additional states component for the tip states)
#' @param observed_tipstates a vector of tip states where each element corresponds to a tip; the name attribute contains the tip labels
#' @param states a vector of states of the discrete character
#' @return a tree with two additional components, node.states_parsimony_observed and node.states_parsimony_simulated, attached
sim_par <- function(tree, observed_tipstates, states) { # only work for a single site/character
  
  # parsimony history for the observed tip states
  pr_oberserved <- phangorn::ancestral.pars(tree = tree, 
                                            data = phyDat(t(t(observed_tipstates)), type = "USER", levels = states), 
                                            type = "MPR", return = "phyDat")
  
  # parsimony history for the simulated tip states
  simulated_tipstates <- tree$states
  simulated_tipstates[names(observed_tipstates)[observed_tipstates == "?"]] <- "?"
  pr_simulated <- phangorn::ancestral.pars(tree = tree, 
                                           data = phyDat(t(t(simulated_tipstates)), type = "USER", levels = states), type = "MPR", return = "phyDat")
  
  tree$node.states_parsimony_observed <- matrix(states[as.integer(pr_oberserved[tree$edge])], ncol = 2)
  tree$node.states_parsimony_simulated <- matrix(states[as.integer(pr_simulated[tree$edge])], ncol = 2)
  
  return(tree)
}

#' Reconstruct some number of parsimony histories of a discrete-character over a list of bifurcation trees conditioning on the observed tip states and the simulated tip states, respectively
#' @param trees a list of bifuracting trees (each of them is a phylo object with an additional states component for the tip states)
#' @param observed_tipstates a vector of tip states where each element corresponds to a tip; the name attribute contains the tip labels
#' @param states a vector of states of the discrete character
#' @param indices generation indices that will only used for printing progress to screen
#' @return a list of trees where each has two additional components, node.states_parsimony_observed and node.states_parsimony_simulated, attached
sim_pars <- function(trees, observed_tipstates, states, indices = NULL) {
  
  ntrees <- length(trees)
  verbose <- !is.null(indices)
  for (l in 1:ntrees) {
    if (verbose && l %% 50 == 0) {
      cat(paste0("parsimony reconstruction no. ", indices[l], ".\n"))
    }
    if (!any(is.na(trees[[l]]))) {
      trees[[l]] <- sim_par(tree = trees[[l]], observed_tipstates, states)
    }
  }
  
  return(trees)
}

#' Reconstruct parsimony histories of a discrete-character over a list of bifurcation trees conditioning on the observed tip states and the simulated tip states, respectively
#' @param trees a list of bifuracting trees (each of them is a phylo object with an additional states component for the tip states)
#' @param observed_tipstates a vector of tip states where each element corresponds to a tip; the name attribute contains the tip labels
#' @param states a vector of states of the discrete character
#' @param ncores number of computer cores to use (if more than one then simulations will be parallelized)
#' @return a list of trees where each has two additional components, node.states_parsimony_observed and node.states_parsimony_simulated, attached
simulate_pars <- function(trees, observed_tipstates, states, ncores = 1L) {
  
  ntrees <- length(trees)
  if (ncores > 1L && ntrees > 50L) { # multi-core parallelization
    
    nchunks <- ncores * 4L
    if (ntrees < nchunks) {
      nchunks <- ntrees
    }
    nhis_percore <- floor(ntrees / nchunks)
    his_indices <- vector("list", nchunks)
    k <- 0L
    for (j in 1:nchunks) {
      if (j < nchunks) {
        his_indices[[j]] <- 1L:nhis_percore + k
      } else {
        his_indices[[j]] <- (k + 1L):ntrees
      }
      k <- k + nhis_percore
    }
    
    trees_all <- pblapply(his_indices, function(his_idx) sim_pars(trees = trees[his_idx], observed_tipstates, states, indices = his_idx), cl = ncores)
    trees <- do.call(c, trees_all)
  } else { # single core
    trees <- sim_pars(trees, observed_tipstates, states, indices = seq_along(trees))
  }
  
  return(trees)
}

#' Compute the piecewise constant parsimony statistics
#' @param tree a bifuracting tree of class phylo
#' @param epoch_bound_ages a vector containing the age of time interval boundaries for a piecewise constant geographic model 
#' (ordered decreasingly so the first two elements are the boundaries of the most ancient interval)
#' @param node.states identical to the edge matrix of a phylo object except here each element is the node state instead of node index
#' @return a vector (each element correspond to a time interval, from present to ancient) of piecewise constant parsimony statistics
parsimonyscore_epoch_calculator <- function(tree, epoch_bound_ages, node.states) {

  # make sure the boundaries are orderred decreasingly
  epoch_bound_ages <- sort(epoch_bound_ages, decreasing = T)
  epoch_num <- length(epoch_bound_ages) - 1L
  
  nedges <- nrow(tree$edge)
  node_times <- ape::node.depth.edgelength(tree)
  node_ages <- max(node_times) - node_times
  
  pieceduration_mat <- matrix(0, ncol = epoch_num, nrow = nedges, byrow = T)
  for (k in 1:nedges) { # loop over branches to get the interval mapping matrix
    
    node_anc <- tree$edge[k, 1]
    node_dec <- tree$edge[k, 2]
    age_anc <- node_ages[node_anc]
    age_dec <- node_ages[node_dec]
    
    for (m in 1:epoch_num) {
      piece_duration <- min(age_anc, epoch_bound_ages[m]) - max(age_dec, epoch_bound_ages[m + 1])
      if (piece_duration > 0) {
        pieceduration_mat[k, m] <- piece_duration
      }
    }
  }
  
  nchanges_epoch <- integer(epoch_num)
  for (k in 1:nedges) { # loop over branches again to count the number of events
    if (length(unique(node.states[k, ])) > 1 && all(node.states[k, ] != "?")) { # when the start and end states are different (therefore a change must had occurred)
      pieces_duration <- pieceduration_mat[k, ]
      
      if (sum(pieces_duration > 0) > 1) {
        epoch_id <- sample.int(epoch_num, size = 1, prob = pieces_duration)
        nchanges_epoch[epoch_id] <- nchanges_epoch[epoch_id] + 1L
      } else if (sum(pieces_duration > 0) == 1) {
        nchanges_epoch[pieces_duration > 0] <- nchanges_epoch[pieces_duration > 0] + 1L
      } else {
        stop("branch cannot be zero length when there is a change on it.\n")
      }
    }
  }
  
  return(rev(nchanges_epoch))
}

#' computing a posterior predictive p-value of a posterior-predictive distribution
#' @param simulated distribution (a vector) of simulated values
#' @param observed distribution (a vector) of observed values
#' if only one vector is provided then it would be deemed as a distribution of simulated - observed values
#' @param method method used to compute the posterior-predictive p-value: 
#' the default option ("doubletail") computes the conventional two-tailed p-value (a value that is very close to 0 or to 1 indicates inadequacy),
#' while the other options compute a single-tailed p-value. The "symmetric" option assumes the distribution is symmetric (so that the left tail and right tail are identical).
#' The last two options compute a probability density from the vector of simulated - observed values first, 
#' and then calculate the fraction of these values with a probability density that is smaller than the probability density of zero (under option "montecarlo"), or
#' calculate the area under the probability density curve between the intervals whose probability density is smaller than the probability density of zero (under option "integrate").
#' @return a posterior-predictive p-value
get_ppp <- function(simulated, observed = 0, method = "doubletail") {
  if (!identical(observed, 0)) {
    sim_minus_obs <- simulated - observed
  } else {
    sim_minus_obs <- simulated
  }
  if (min(sim_minus_obs) > 0 || max(sim_minus_obs) < 0) return(0)
  
  if (method == "doubletail") {
    ppp <- (sum(sim_minus_obs > 0) + sum(sim_minus_obs == 0) / 2) / length(sim_minus_obs)
    return(ppp)
  }
  
  if (method == "symmetric") {
    sim_minus_obs_ecdf <- ecdf(sim_minus_obs)
    ppp <- sim_minus_obs_ecdf(0)
    if (ppp > 0.5) ppp <- 1 - ppp
    return(ppp * 2)
  }
  
  nsims <- length(sim_minus_obs)
  sim_minus_obs_tab <- table(sim_minus_obs)
  sim_minus_obs_unique <- as.integer(names(sim_minus_obs_tab))
  
  sim_minus_obs_dens <- density(sim_minus_obs, bw = "SJ", from = min(sim_minus_obs_unique), to = max(sim_minus_obs_unique))
  sim_minus_obs_pdf <- approxfun(x = sim_minus_obs_dens$x, y = sim_minus_obs_dens$y)
  
  sim_minus_obs_unique_pd <- sim_minus_obs_pdf(sim_minus_obs_unique)
  zero_pd <- sim_minus_obs_pdf(0)
  
  if (method == "montecarlo") {
    nsims_moreextreme <- sum(sim_minus_obs_tab[sim_minus_obs_unique_pd < zero_pd])
    ppp <- nsims_moreextreme / nsims
  } else if (method == "integrate") {
    moreextreme_area <- integrate(function(x) ifelse(sim_minus_obs_pdf(x) < zero_pd, sim_minus_obs_pdf(x), 0), 
                                  lower = min(sim_minus_obs_unique), upper = max(sim_minus_obs_unique), subdivisions = 1e5,
                                  abs.tol = 0.01)
    ppp <- moreextreme_area$value
  }
  
  return(ppp)
}
