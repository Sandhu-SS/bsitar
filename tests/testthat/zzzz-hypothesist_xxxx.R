

# full.argsx

# marginaleffects_hypothesis_args[['by']] <- full.argsx[['by']]
# marginaleffects_hypothesis_args[['hypothesis']] <- difference ~ pairwise

# devtools::load_all()




full.argsx2 <- full.argsx

# hypothesis_args[['by']] <- 'sex'
 full.argsx2[['hypothesis']] <- difference ~ pairwise | sex

dtttt <- get_comparison_hypothesis(draws_list_dtx, full.argsx2,by=c( 'sex'), 
                                   evaluate_comparison = T,
                                   evaluate_hypothesis = T,
                                  get_eqpd_form = F, format = T,
                                 get_eqpd = F,
                                 )
