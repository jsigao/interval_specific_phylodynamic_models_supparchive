# this script contains functions that intake sampled history or dated phylogeny and then loop over time intervals to
# count the number of events or number of active viral lineages in each arbitrarily specified time interval

#' Extract and store the pairwise number of events from a BEAST history log-file output
#' @param log_path path to the BEAST history log file
#' @param burnin fraction of or absolute number of generations of the log file to be discarded as burnin
#' @return a list of matrix (each element corresponds to a pair) of list (each element corresponds to a generation) 
#' of data frames (each row corresponds to an event)
get_complete_history <- function(log_path, burnin = 0) {
  
  # get the number of states
  xml_path <- grep("MLE", list.files(dirname(log_path), recursive = T, full.names = T, pattern = "*.xml$"), invert = T, value = T)[1]
  x <- scan(file = xml_path, what = character(), sep = "\n", strip.white = F, blank.lines.skip = F)
  state.names <- gsub("^\\s+|\\s+$", "", gsub("<state code|=|\t|\"|/>", "", x[grep("<state code", x)]))
  state.num <- length(state.names)
  
  history_log <- read.table(log_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  if (burnin > 0) {
    if (burnin < 1) {
      burnin <- ceiling(burnin * nrow(history_log))
    }
    history_log <- history_log[-(1:burnin), ]
  }
  
  # put the transitions in each generation into a data frame
  history_df <- vector("list", nrow(history_log))
  for (l in 1:length(history_df)) {
    history <- unlist(strsplit(history_log$completeHistory_1[l], "\\{\\{|\\}\\,\\{|\\}\\}"))
    history <- history[-c(1, length(history))]
    history <- data.frame(matrix(unlist(strsplit(history, ",")), ncol = 4, byrow = T), stringsAsFactors = F)
    history <- history[, -1]
    
    colnames(history) <- c("age", "from", "to")
    history$age <- as.numeric(history$age)
    history_df[[l]] <- history
  }
  
  rm(history_log)
  gc()
  
  # fetch the events between each pair and then store them in the corresponding cell of the large matrix
  history_eachpair <- matrix(list(), nrow = state.num, ncol = state.num, byrow = T)
  for (i in 1:state.num) {
    for (j in 1:state.num) {
      if (i != j) {
        
        history_thispair <- vector("list", length(history_df))
        for (k in 1:length(history_df)) {
          history_thispair[[k]] <- history_df[[k]][history_df[[k]]$from == state.names[i] & history_df[[k]]$to == state.names[j], ]
        }
        
        history_eachpair[[i, j]] <- history_thispair
      }
    }
  }
  rownames(history_eachpair) <- state.names
  colnames(history_eachpair) <- state.names
  
  rm(history_df)
  gc()
  
  return(list(history_eachpair = history_eachpair))
}


#' Count the number of events between each pair of areas in each time interval
#' @param history_eachpair a matrix of lists where each list corresponds to a given pair of areas, 
#' and the element of each list corresponds to a history sample in the posterior 
#' (specifically, it is a data frame containing the events occurred in the given history)
#' @param epoch_bounds boundaries (a vector of numbers) delineating the time intervals; e.g., if the entire sampling period is 
#' 100 days old and we want to count the number of events per day, then this vector needs to be 0:100
#' @return a list of matrices, where each matrix corresponds to a time interval (the first matrix represents the last time interval); 
#' the element of each matrix is a vector of number of events occurred between the associated pair of areas in that time interval,
#' and the element of this vector corresponds to that statistic in each sampled history
get_nevents_perinterval <- function(history_eachpair, epoch_bounds) {
  
  epoch_num <- length(epoch_bounds)
  state_name <- rownames(history_eachpair)
  state_num <- ncol(history_eachpair)
  
  counts_epoch <- vector("list", epoch_num)
  for (m in 1:epoch_num) {
    
    counts <- matrix(list(), nrow = state_num, ncol = state_num, byrow = T)
    
    for (k in 1:state_num) {
      for (l in 1:state_num) {
        if (k != l) {
          
          count <- integer(length(history_eachpair[[k, l]]))
          for (n in 1:length(history_eachpair[[k, l]])) {
            
            if (nrow(history_eachpair[[k, l]][[n]]) > 0) {
              if (m < epoch_num) {
                count[n] <- sum(history_eachpair[[k, l]][[n]]$age < epoch_bounds[m + 1] & history_eachpair[[k, l]][[n]]$age >= epoch_bounds[m])
              } else if (m == epoch_num) {
                count[n] <- sum(history_eachpair[[k, l]][[n]]$age >= epoch_bounds[m])
              }
            }
          }
          
          counts[[k, l]] <- count
        }
      }
    }
    
    rownames(counts) <- state_name
    colnames(counts) <- state_name
    counts_epoch[[m]] <- counts
  }
  
  return(counts_epoch)
}


#' Count the number of events occurred on a given dispersal route in each time interval
#' @param counts_epoch a list of matrices, where each matrix corresponds to a time interval (the first matrix represents the last time interval); 
#' the element of each matrix is a vector of number of events occurred between the associated pair of areas in that time interval,
#' and the element of this vector corresponds to that statistic in each sampled history
#' @param origin the origin area(s) of the dispersal route (can be either an area index or a vector of area indices)
#' @param dest the destination area(s) of the dispersal route (can be either an area index or a vector of area indices)
#' @return a list of vectors, where each vector corresponds to a time interval (the first vector represents the last time interval);
#' each element of each vector is the number of events occurred on a given dispersal route in the given time interval in the corresponding sample of history
get_nevents_oneroute_perinterval <- function(counts_epoch, origin, dest) {
  
  epoch_num <- length(counts_epoch)
  state_num <- ncol(counts_epoch[[1]])
  state_name <- colnames(counts_epoch[[1]])
  
  nhis <- length(counts_epoch[[1]][[1, 2]])
  nevents_epoch <- vector("list", epoch_num)
  
  for (m in 1:epoch_num) {
    nevents <- integer(nhis)
    for (n in 1:nhis) {
      nevents[n] <- sum(sapply(counts_epoch[m], 
                               function(x) sum(unlist(sapply(x[row(x) %in% origin & col(x) %in% dest], "[[", n)))))
    }
    nevents_epoch[[m]] <- nevents
  }
  
  return(nevents_epoch)
}
