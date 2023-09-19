

get_idata <- function(newdata, idVar, timeVar, times = NULL, length.out = 10, xrange = 1) {
  all_times <- TRUE
  if(is.null(xrange)) xrange <- 1 else xrange <- xrange
  times_orig <- newdata[[timeVar]]
  times_orig <- times_orig[!is.na(times_orig)]
  
  if (is.null(times) || !is.numeric(times)) {
    times <- seq(min(times_orig), max(times_orig), length.out = length.out)
  }
  
  # this is when no random effects and this groupvar is NULL
  # therefore, an artificial group var created 
  # check utils-helper function lines 60
  
  if(nlevels(newdata[[idVar]]) == 1) {
    newdata <- newdata %>% distinct(newdata[[timeVar]], .keep_all = T) %>% 
      arrange(!!as.name(timeVar))
  }
  
  
  id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
  
  # , setminmax_xrange = NULL
  # if(is.null(setminmax_xrange)) {
  #   last_time <- tapply(newdata[[timeVar]], id, max)
  #   first_time <- tapply(newdata[[timeVar]], id, min)
  # }
  
  last_time <- tapply(newdata[[timeVar]], id, max)
  first_time <- tapply(newdata[[timeVar]], id, min)
  
  newdata_nomiss <- newdata[complete.cases(newdata), ]
  id_nomiss <- match(newdata_nomiss[[idVar]], unique(newdata_nomiss[[idVar]]))
  n <- length(unique(id_nomiss))
  
  if(xrange == 1) {
    times_to_pred <- list()
    for (i in 1:length(unique(newdata[[idVar]]))) {
      numx <- as.character(i)
      times_to_pred[[numx]] <- seq(first_time[i], last_time[i], length.out = length.out)
    }
  }
  
  if(xrange == 2) {
    times_to_pred <- lapply(last_time, function (t) 
      if (all_times) times else times[times > t])
  }
  
  id_pred <- rep(seq_len(n), sapply(times_to_pred, length))
  
  right_rows <- function (data, times, ids, Q_points) {
    fids <- factor(ids, levels = unique(ids))
    if (!is.list(Q_points))
      Q_points <- split(Q_points, row(Q_points))
    ind <- mapply(findInterval, Q_points, split(times, fids))
    ind[ind < 1] <- 1
    rownams_id <- split(row.names(data), fids)
    ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
    data[c(ind), ]
  }
  
  newdata_pred <- right_rows(newdata, newdata[[timeVar]], id, times_to_pred)
  
  newdata_pred[[timeVar]] <- unlist(times_to_pred)
  
  newdata_pred
}



# v = seq(6, 20, by = .10)

# get_idata %>% newinterpolate(., idVar = 'id', timeVar = 'age', times = NULL, xrange = 2)