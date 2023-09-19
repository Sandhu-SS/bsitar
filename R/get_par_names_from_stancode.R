

get_par_names_from_stancode <- function(code, 
                                        full = TRUE, 
                                        section =  'parameters',
                                        what = '') {
  regex_for_section <- paste(".*(",section,"\\s*\\{.*?\\}).*", sep = '')
  filtered_stan_code <- gsub(code, pattern = regex_for_section, replacement = "\\1")
  
  zz <- strsplit(filtered_stan_code, "\n")[[1]][-1]
  collect <- c()
  collect_full <- c()
  for (i in 1:length(zz)-1) {
    if(!(identical(zz[i], character(0))))  {
      # print(zz[i])
      t <- sub(";.*", "", zz[i])
      t_full <- t
    #  t_full <- gsub("[[:space:]]", "", t_full)
      t_full <- gsub("^ *|(?<= ) | *$", "", t_full, perl=T)
      if(what == "") {
        get_t_full <- t_full
        collect_full <- c(collect_full, get_t_full)
      } else if(what != "") {
        get_t_full <- t_full
        get_t_full <- get_t_full[grepl(paste0(what, "_"), get_t_full)]
        collect_full <- c(collect_full, get_t_full)
      }
      t <- tail(strsplit(t,split=" ")[[1]],1)
      collect <- c(collect, t)
    }
  }
  if(!full) {
    out <- collect
  }
  if(full) {
    out <- collect_full
  }
  out
}


# xscode <- stancode(female_1444)
# get_par_names_from_stancode(xscode, full = T, what = '') # sd L z
