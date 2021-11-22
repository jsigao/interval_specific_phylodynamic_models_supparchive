#' sort a vector with ties broken randomly
#' @param x a vector to be sorted
#' @param decreasing whether sorting the vector decreasingly or increasingly
#' @return a sorted vector
sort_tierandom <- function(x, decreasing = F) {
  if (decreasing) {
    x[match(seq_along(x), rank(-x, ties.method = "random"))]
  } else {
    x[match(seq_along(x), rank(x, ties.method = "random"))]
  }
}

#' calculate a Highest Probability Density (HPD) Interval of a distribution
#' @param x a vector of samples
#' @param prob size of the HPD (a value between zero and 1)
#' @param method either treating the samples as integers (method = "integer") or discrete variables 
#' (method = "discrete"; i.e., if samples are integer then 1 is equally close to 2 and 3). 
#' @return an HPD set computed from the provided vector samples corresponding to the size specified by prob
hpd <- function(x, prob, method = "discrete") {
  
  if (method == "discrete") {
    freq <- table(x) / length(x)
    freq <- sort_tierandom(freq, decreasing = T)
    x_unique <- as.integer(names(freq))
  } else if (method == "integer") {
    x_unique <- unique(x)
    bw <- bw.nrd0(x)
    # if (bw < 0.1) bw <- 0.1
    x_dens <- density(x, bw = bw, from = min(x_unique), to = max(x_unique))
    x_pdf <- approxfun(x = x_dens$x, y = x_dens$y)
    
    x_unique_inter <- min(x_unique):max(x_unique)
    freq_inter <- x_pdf(x_unique_inter)
    freq_inter <- freq_inter / sum(freq_inter)
    names(freq_inter) <- x_unique_inter
    freq_inter <- sort_tierandom(freq_inter, decreasing = T)
    x_unique <- as.integer(names(freq_inter))
    freq <- freq_inter
  } else {
    stop("the specified hpd method has not been implemented")
  }

  nvalues <- length(freq)
  prob_left <- prob
  k <- 0L
  for (i in 1:nvalues) {
    prob_inc <- prob_left / freq[i]
    if (prob_inc >= 1 || runif(1L) <= prob_inc) {
      k <- k + 1L
    }
    prob_left <- prob_left - freq[i]
    if (prob_left <= 0) break
  }
  
  if (k == 0) return(integer())
  res <- x_unique[seq_len(k)]
  return(res)
}

#' obtain the number of dispersal events from a simulated history list (from the true values)
#' @param histories_true a list of simulated histories where each one is a phylo and simmap object that contains both the phylogeny and the dispersal history
#' @param nsims the first nsims number of simulated histories to include
#' @return a list where each element (corresponding to one simulated dataset) is in turn a vector (where each matrix corresponds to a sample of MCMC) that contains pairwise number of dispersal events
get_counts_true <- function(histories_true, nsims = 0) {
  if (nsims == 0 || nsims > length(histories_true)) {
    nsims <- length(histories_true)
  }
  
  states <- sort(unique(histories_true[[1]]$states))
  nstates <- length(states)
  nroutes <- choose(nstates, 2) * 2
  nparams <- 1L + nroutes
  
  nevents_true_all <- sapply(histories_true, function(x) sum(lengths(x$maps)) - length(x$maps))
  params_true_all <- vector("list", nsims)
  for (j in 1:nsims) {
    params_true_all[[j]] <- integer(nparams)
    params_true_all[[j]][1] <- nevents_true_all[j]
  }
  
  for (j in 1:nsims) {
    maps <- histories_true[[j]]$maps
    maps <- maps[lengths(maps) > 1]
    maps <- lapply(maps, names)
    maps_str <- sapply(maps, function(x) paste0("+", paste(x, collapse = "+"), "+"))
    m <- 2L
    
    if (length(maps_str) > 0) {
      for (k in 1:nstates) {
        for (l in 1:nstates) {
          if (k != l && length(maps_str) > 0) {
            state_str <- paste(states[c(k, l)], collapse = "+")
            state_strpp <- paste0("+", state_str, "+")
            params_true_all[[j]][m] <- sum(maps_str == state_strpp)
            maps_str <- maps_str[maps_str != state_strpp]
            maps_str_this <- grep(state_strpp, maps_str, fixed = T, value = T)
            
            if (length(maps_str_this) > 0) {
              params_true_all[[j]][m] <- params_true_all[[j]][m] + sum(lengths(strsplit(maps_str_this, state_str, fixed = T)) - 1)
            }
            
            m <- m + 1L
          }
        }
      }
    }
  }
  
  return(params_true_all)
}

#' obtain the distribution of the number of dispersal events from a parameter list (for the estimated values from simulated datasets)
#' @param params_allsims a list where each element (corresponding to one simulated dataset) is a list returned by get_BF_counts function that contains the samples from a BEAST log file
#' @return a list where each element (corresponding to one simulated dataset) is in turn a list of matrices (where each matrix corresponds to a sample of MCMC) that contains pairwise number of dispersal events
get_counts_sim <- function(params_allsims) {
  nsims <- length(params_allsims)
  counts_mat_all <- vector("list", nsims)
  nepochs <- ifelse("counts_mat_all" %in% names(params_allsims[[1]][[1]]), length(params_allsims[[1]]), length(params_allsims[[1]][[1]]))
  
  for (j in 1:nsims) {
    if (nepochs == 1) {
      counts_mat_all[[j]] <- params_allsims[[j]][[1]]$counts_mat_all
    } else {
      counts_mat_all[[j]] <- params_allsims[[j]][[1]][[1]]$counts_mat_all
      nsamples <- length(params_allsims[[j]][[1]][[1]]$counts_mat_all)
      for (k in 2:nepochs) {
        for (l in 1:nsamples) {
          counts_mat_all[[j]][[l]] <- counts_mat_all[[j]][[l]] + params_allsims[[j]][[1]][[k]]$counts_mat_all[[l]]
        }
      }
    }
  }
  
  return(counts_mat_all)
}

#' compute the HPD for a distribution of the number of dispersal events 
#' @param counts_sim a list where each element (corresponding to one simulated dataset) is in turn a list of matrices (where each matrix corresponds to a sample of MCMC) that contains pairwise number of dispersal events
#' @param probs vector of sizes of the HPD (each is a value between zero and 1)
#' @param method either treating the samples as integers (method = "integer") or discrete variables 
#' (method = "discrete"; i.e., if samples are integer then 1 is equally close to 2 and 3).
#' @return a list where each of its elements corresponds to the number of dispersal events among all areas (the first element) or over a given dispersal route (the remaining elements);
#' each such element is in turn a list where each element corresponds to a simulated dataset; this list then contains HPD sets of the number of dispersal events
#' where each set corresponds to a given probability provided in the probs vector
get_counts_sim_hpd <- function(counts_sim, probs, hpd_method = "discrete") {
  nstates <- nrow(counts_sim[[1]][[1]])
  nroutes <- choose(nstates, 2) * 2
  nparams <- 1L + nroutes
  nsims <- length(counts_sim)
  
  params_credsets <- vector("list", nparams)
  for (i in 1:nparams) {
    params_credsets[[i]] <- vector("list", nsims)
  }
  
  for (j in 1:nsims) {
    counts_total <- sapply(counts_sim[[j]], sum)
    params_credsets[[1]][[j]] <- lapply(probs, function(prob) hpd(counts_total, prob, method = hpd_method))
  }
  
  m <- 2L
  for (k in 1:nstates) {
    for (l in 1:nstates) {
      if (k != l) {
        for (j in 1:nsims) {
          counts_kl <- sapply(counts_sim[[j]], function(y) y[k, l])
          params_credsets[[m]][[j]] <- lapply(probs, function(prob) hpd(counts_kl, prob, method = hpd_method))
        }
        m <- m + 1L
      }
    }
  }
  
  return(params_credsets)
}

#' compute mean of a distribution of the number of dispersal events
#' @param counts_sim a list where each element (corresponding to one simulated dataset) is in turn a list of matrices (where each matrix corresponds to a sample of MCMC) that contains pairwise number of dispersal events
#' @return a list of vectors where each vector contains the mean estimate of the number of dispersal events among all areas (the first vector) 
#' or over a given dispersal route (the remaining vectors) inferred from each simulated dataset (i.e., each element of a vector corresponds to a simulated dataset)
get_counts_sim_mean <- function(counts_sim) {
  nstates <- nrow(counts_sim[[1]][[1]])
  nroutes <- choose(nstates, 2) * 2
  nparams <- 1L + nroutes
  
  params_mean <- vector("list", nparams)
  params_mean[[1]] <- sapply(counts_sim, function(x) mean(sapply(x, sum)))
  
  m <- 2L
  for (k in 1:nstates) {
    for (l in 1:nstates) {
      if (k != l) {
        params_mean[[m]] <- sapply(counts_sim, function(x) mean(sapply(x, function(y) y[k, l])))
        m <- m + 1L
      }
    }
  }
  
  return(params_mean)
}

#' obtain the distribution of a given parameter from a parameter list (for the estimated values from simulated datasets) generated by get_BF_counts
#' @param param is a list returned by get_BF_counts function that contains the samples from a BEAST log file
#' @param param_name name of the parameter (as a member of the param list) to summarize
#' @return a list where each element (corresponding to one time interval) is in turn a list where each element (corresponding to a sample of MCMC) is a vector of estimated parameter values
get_param <- function(params, param_name) {
  if (param_name == "mu_all") {
    nQepochs <- ifelse(param_name %in% names(params[[1]]), length(params), length(params[[1]]))
    if (nQepochs == 1) {
      nmuepochs <- ifelse(is.list(params[[1]][[param_name]]), length(params[[1]][[param_name]]), 1L)
    } else {
      nmuepochs <- 0
      for (k in 1:nQepochs) {
        nmuepochs <- nmuepochs + ifelse(is.list(params[[1]][[k]][[param_name]]), length(params[[1]][[k]][[param_name]]), 1L)
      }
    }
    
    param <- vector("list", nmuepochs)
    if (nmuepochs == 1) {
      param[[1]] <- params[[1]][[param_name]]
    } else {
      if (nQepochs == 1) {
        param <- params[[1]][[param_name]]
      } else {
        j <- 1L
        for (k in 1:nQepochs) {
          nmuepochs_this <- 1L
          if (is.list(params[[1]][[k]][[param_name]])) {
            nmuepochs_this <- length(params[[1]][[k]][[param_name]])
            param[j - 1 + 1:nmuepochs_this] <- params[[1]][[k]][[param_name]]
          } else {
            param[[j]] <- params[[1]][[k]][[param_name]]
          }
          
          if (k > 1) {
            if (identical(param[[j]], param[[j - 1]])) {
              param[[j]] <- NULL
              j <- j - 1L
            }
          }
          j <- j + nmuepochs_this
        }
        
        param <- param[1:(j - 1)]
      }
    }
    
  } else {
    nQepochs <- ifelse(param_name %in% names(params[[1]]), length(params), length(params[[1]]))
    param <- vector("list", nQepochs)
    if (nQepochs == 1) {
      param[[1]] <- params[[1]][[param_name]]
    } else {
      for (k in 1:nQepochs) {
        param[[k]] <- params[[1]][[k]][[param_name]]
      }
    }
  }
  
  return(param)
}
