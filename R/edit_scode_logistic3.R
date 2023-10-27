
#' An internal function to edit stancode for tripple logistic model
#' @param stancode A string character of stan code
#' @param genq_only A logical (default \code{FALSE}) to indicate whether to 
#' return only the generated quantity sub code.
#' @param normalize A logical (default \code{TRUE}) to indicate whether to 
#' include the normalizing constant in the prior target density.
#' @keywords internal
#' @return A character string.
#' @noRd
#'
edit_scode_logistic3 <- function(stancode, 
                                 genq_only = FALSE, 
                                 normalize = TRUE) {
  
  # Rename transformed parameters and parameters for ease of processing
  
  true_name_p   <- 'parameters'
  true_name_tp  <- 'transformed parameters'
  true_name_td  <- 'transformed data'
  
  tempt_name_p  <- 'xxxxxxxxxxx'
  tempt_name_tp <- 'zzzzzzzzzzzz'
  tempt_name_td <- 'ddddddddddd'
  
  clines_tp <- get_par_names_from_stancode(stancode,
                                           section =  true_name_tp,
                                           semicolan = TRUE,
                                           full = TRUE)
  
  clines_p <- get_par_names_from_stancode(stancode,
                                          section =  true_name_p,
                                          semicolan = TRUE,
                                          full = TRUE)
  
  clines_m <- get_par_names_from_stancode(stancode,
                                          section =  'model',
                                          semicolan = TRUE,
                                          full = TRUE)
  
  
  editedcode    <- stancode 
  editedcode    <- gsub(true_name_tp, tempt_name_tp, editedcode, fixed = T)
  editedcode    <- gsub(true_name_p,  tempt_name_p,  editedcode, fixed = T)
  editedcode    <- gsub(true_name_td,  tempt_name_td,  editedcode, fixed = T)
  
  editedcode2 <- editedcode
  
  
  clines_tp2 <- c()
  for (il in clines_tp) {
    il <- gsub(pattern = "//", replacement = "//", x = il, fixed = T)
    il <- gsub(pattern = "//[^\\\n]*", replacement = "", x = il)
    if(!grepl('^lprior', gsub_space(il)) & 
       !grepl('^reallprior', gsub_space(il)) & 
       !grepl('^-', gsub_space(il))) {
      if(!is_emptyx(il)) {
        clines_tp2 <- c(clines_tp2, il) 
      }
    }
  }
  
  clines_tp <- clines_tp2
  
  
  clines_p2 <- c()
  for (il in clines_p) {
    il <- gsub(pattern = "//", replacement = "//", x = il, fixed = T)
    il <- gsub(pattern = "//[^\\\n]*", replacement = "", x = il)
    if(!grepl('^lprior', gsub_space(il)) & 
       !grepl('^reallprior', gsub_space(il)) & 
       !grepl('^-', gsub_space(il))) {
      if(!is_emptyx(il)) {
        clines_p2 <- c(clines_p2, il) 
      }
    }
  }
  
  clines_p <- clines_p2
  
  
  
  b_what_by_pair <- matrix(NA, 9, 2)
  K_name <- "K"
  move_to_tp <- add_move_to_tp <- c()
  for (clines_tpi in clines_p) {
    parameter_name_temp <- sub('.+](.+)', '\\1', clines_tpi)
    parameter_name_temp <- gsub(";", "", parameter_name_temp)
    parameter_name_temp <- gsub_space(parameter_name_temp)
    for (igr in 1:9) {
      if(grepl(paste0("^b", "_"), parameter_name_temp) & 
         grepl(paste0("_", letters[igr]), parameter_name_temp)) {
        move_to_tp <- c(move_to_tp, clines_tpi)
        if(letters[igr] == 'a')  poparm <- paste0("raw_d[", 1, "]")
        if(letters[igr] == 'b')  poparm <- paste0("raw_v[", 3, "]")
        if(letters[igr] == 'c')  poparm <- paste0("raw_t[", 1, "]")
        if(letters[igr] == 'd')  poparm <- paste0("raw_d[", 3, "]")
        if(letters[igr] == 'e')  poparm <- paste0("raw_v[", 1, "]")
        if(letters[igr] == 'f')  poparm <- paste0("raw_t[", 2, "]")
        if(letters[igr] == 'g')  poparm <- paste0("raw_d[", 2, "]")
        if(letters[igr] == 'h')  poparm <- paste0("raw_v[", 2, "]")
        if(letters[igr] == 'i')  poparm <- paste0("raw_t[", 3, "]")
        
        add_move_to_tp_k_dim <- paste0(K_name, "_", letters[igr])
        add_move_to_tp_ <- paste0(" = ", poparm, " + ", "rep_vector(0.0, ",
                                  add_move_to_tp_k_dim, ")")
        add_move_to_tp_2 <- gsub(";", paste0(add_move_to_tp_, ";"),  clines_tpi, fixed = T)
        
        add_move_to_tp <- c(add_move_to_tp, add_move_to_tp_2)
        parameter_name <- sub('.+](.+)', '\\1', clines_tpi)
        parameter_name <- gsub(";", "", parameter_name)
        parameter_name <- gsub_space(parameter_name)
        b_what_by_pair[igr , 1] <- parameter_name
        b_what_by_pair[igr , 2] <- poparm
      } 
    }
  } 
  # b_what_by_pairx <<- b_what_by_pair
  
  keep_clines_p <- setdiff(clines_p, move_to_tp)
  
  add_keep_clines_p <- "
  simplex[3+1] simp_d;
  simplex[3+1] simp_v;
  simplex[3+1] simp_t;
  "
  
  keep_clines_p <- c(add_keep_clines_p, keep_clines_p)
  prepare_p <- paste(keep_clines_p, collapse = "\n") 
  
  add_keep_clines_tp <- 
    "positive_ordered [3] raw_d;
  positive_ordered [3] raw_v;
  positive_ordered [3] raw_t;
  raw_d = min_d + head(cumulative_sum(simp_d), 3) * (max_d - min_d);
  raw_v = min_v + head(cumulative_sum(simp_v), 3) * (max_v - min_v);
  raw_t = min_t + head(cumulative_sum(simp_t), 3) * (max_t - min_t);"
  
  
  move_to_tp <- c(add_keep_clines_tp, add_move_to_tp)
  move_to_tp <- paste(move_to_tp, collapse = "\n") 
  
  
  move_to_td <- 
    "real min_d = 0;
  real max_d = 200; 
  real min_v = 0;
  real max_v = 5; 
  real min_t = 0;
  real max_t = 18;"
  
  move_to_td <- paste(move_to_td, collapse = "\n") 
  
  # cat(prepare_p)
  # cat(move_to_tp)
  # print(b_what_by_pair)
  # stop()
  
  
  
  
  # Add space to model block elements
  # zz_c <- c()
  # for (iz in move_to_m) {
  #   zz_c <- c(zz_c, paste0("  ", iz))
  # }
  # move_to_m <- paste(zz_c, collapse = '\n')
  
  
  
  
  
  
  for (il in clines_p) {
    editedcode2 <- gsub(pattern = "//", replacement = "//", x = editedcode2, fixed = T)
    editedcode2 <- gsub(pattern = "//[^\\\n]*", replacement = "", x = editedcode2)
    editedcode2 <- gsub(paste0(il, ""), "", editedcode2, fixed = T)
    
  }
  
  # Remove empty lines
  zz <- strsplit(editedcode2, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    if(!is_emptyx(gsub_space(zz[iz]))) {
      zz_in <- zz[iz]
      zz_c <- c(zz_c, zz_in)
    }
  }
  editedcode2 <- paste(zz_c, collapse = '\n')
  
  
  p_block_syb_by <- paste0("", tempt_name_p, " {")
  p_block_syb_it <- paste0(p_block_syb_by, "\n", prepare_p)
  editedcode2 <- gsub(paste0("", p_block_syb_by), p_block_syb_it, editedcode2, fixed=T, perl=F)
  
  
  p_block_syb_by <- paste0("", tempt_name_tp, " {")
  p_block_syb_it <- paste0(p_block_syb_by, "\n", move_to_tp)
  editedcode2 <- gsub(paste0("", p_block_syb_by), p_block_syb_it, editedcode2, fixed=T, perl=F)
  
  
  p_block_syb_by <- paste0("", tempt_name_td, " {")
  p_block_syb_it <- paste0(p_block_syb_by, "\n", move_to_td)
  editedcode2 <- gsub(paste0("", p_block_syb_by), p_block_syb_it, editedcode2, fixed=T, perl=F)
  
  
  
  
  editedcode2 <- gsub(tempt_name_tp, true_name_tp, editedcode2, fixed = T)
  editedcode2 <- gsub(tempt_name_p,  true_name_p,  editedcode2, fixed = T)
  editedcode2 <- gsub(tempt_name_td,  true_name_td,  editedcode2, fixed = T)
  
  # Replace parameters in priors
  for (igr in 1:nrow(b_what_by_pair)) {
    parameter_name <- b_what_by_pair[igr , 1]
    poparm         <- b_what_by_pair[igr , 2]
    editedcode2 <- gsub(paste0("(", parameter_name, "[", 1, "]"),
                        paste0("(", poparm, ""),
                        editedcode2, fixed = T)
  }
  
  # Remove empty lines
  zz <- strsplit(editedcode2, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    if(!is_emptyx(gsub_space(zz[iz]))) {
      zz_in <- zz[iz]
      zz_c <- c(zz_c, zz_in)
    }
  }
  
  editedcode2 <- paste(zz_c, collapse = '\n')
  
  # cat(editedcode2)
  # stop()
  
  return(editedcode2)
}

