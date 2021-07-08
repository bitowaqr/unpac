
##----------------------------------------------------------
## Multi-step utility aggregation procedure
## Application using EQ-5D 3L data from the 1993 MVH Study 
##  
## Paul Schneider
##  University of Sheffield
##  p.schneider@sheffield.ac.uk
##  27.04.2020
##----------------------------------------------------------


## SETUP -----------------------------------------
  #  rm(list=ls())   # clear global enviroment
  
  # # decimal options
  options(scipen = 99)
  
  # load neccessary packages 
  library(reshape2)
  library(ggplot2)
  library(cowplot)
  
# source convenient functions
  source("./R//support.R")

## READ AND PREPARE THE (TRIMMED) UK EQ-5D 3L MVH DATA SET FOR ANALYSIS-----
  dat = read.csv("./data/mvh_trimmed.csv") # 2,997 MVH participant
  # head(dat)   # see first 5 rows
  dat$healthstate = as.character(dat$healthstate)
  
  # Add 11111 as a health state with a value of 1 for each participant
  dat = dat %>%
    dplyr::add_row(ID=unique(dat$ID),healthstate="11111",tto=1) %>%
    dplyr::arrange(ID)  # sort by ID
  
  # create coefs for regression
  dat = dat %>% bind_cols(hs_to_coef(hs=dat$healthstate))


######## ZERO-DEVISION EXCLUSION 
  # to avoid division by zero, participants with a zero utility range have to be excluded
  non_zero_range = aggregate( tto ~ ID, dat, function(x){ (max(x) - min(x)) > 0 } )
  dat = filter(dat,ID %in% non_zero_range$ID[non_zero_range$tto])
  

### CREATE AN EMPTY SOCIAL VALUE SET: A DATA FRAME WITH ALL 243 HEALTH STATES
  all_states_str = expand.grid(MO=1:3, SC=1:3,UA=1:3, PD=1:3, AD=1:3)
  all_states_str = apply(all_states_str, 1, function(x){paste(x,collapse="")})
  # head(all_states_str)  # health states as strings
  # new df for predictions
  all_states_df = hs_to_coef(all_states_str) # make new df from it
  # head(all_states_df)  # health states in parameter form



## MODEL FITTING ------------------------------- 
  
# The standard approach: Dolan's main effects model
  fit_mean = lm(1-tto ~ MO2 + MO3 + SC2 + SC3 + UA2 + UA3 + PD2 + PD3 + AD2 + AD3 + N3,dat)
  
# Estimating the average social value set
  val_set_mean = 1 - predict(fit_mean,newdata = all_states_df) 
  val_set_mean[1] = 1 # full health = 1 
  val_set_rel_mean = normaliseVS(val_set_mean)
  val_set_mean_rank <- rank(-val_set_mean)

  
# UNPAc 

  # **Step 1: Individual min-max utility normalisation 
  dat_norm = dat %>%
    group_by(ID) %>%
    mutate(
      util_range = max(tto)-min(tto),
      rel_val = (tto - min(tto)) / (max(tto) - min(tto)) 
    )
  
  # ** Step 2: Mean *relative* social value set 
  fit_norm = lm(1-rel_val ~ MO2 + MO3 + SC2 + SC3 + UA2 + UA3 + PD2 + PD3 + AD2 + AD3 + N3, dat_norm)
  # summary(fit_norm)
  
  # ** Step 3a: post-hoc mean anchoring 
  mean_anchoring = mean(dat_norm$util_range)
  
  # ** Step 3b: Re-anchoring 
  fit_norm_mean     = cbind(fit_norm$coefficients,confint(fit_norm)[,1], confint(fit_norm)[,2]) * mean_anchoring
  val_set_norm_mean = predict.anchoredModel(fit_norm_mean, all_states_df)
  val_set_norm_mean[1] = 1
  val_set_rel_norm_mean = normaliseVS(val_set_norm_mean)
  val_set_norm_mean_rank <- rank(-val_set_rel_norm_mean)
  
## COMPARE standard and UNPAc ----
  models <- cbind(coefficients(fit_mean),confint(fit_mean),fit_norm_mean)
  colnames(models) <- c("standard mean","s.l95ci","s.u95ci","UNPAc","u.l95ci","u.u95ci")
  models
  write.csv(models,file = "./output/models.csv",row.names = F)
  
  
  
## COMPARE VALUE SETS -----------
  social_value_sets = cbind(
    "health state" = all_states_str,
    "rank" = val_set_mean_rank,
    "Absolute (relative) value" = paste0(round(val_set_mean,3)," (", round(val_set_rel_mean,3),")"),
    "rank" = val_set_norm_mean_rank,
    "Absolute (relative) value" = paste0(round(val_set_norm_mean,3)," (", round(val_set_rel_norm_mean,3),")"),
    "rank" = formatC(val_set_mean_rank - val_set_norm_mean_rank, digits = 0, format = "f", flag = "+"),
    "Absolute (relative) value" = paste0(
      formatC(val_set_norm_mean - val_set_mean, digits = 2, format = "f", flag = "+"),
      " (", formatC(val_set_rel_norm_mean - val_set_rel_mean, digits = 0, format = "f", flag = "+"),
      ")"
    )
    )
  social_value_sets <- social_value_sets[order(val_set_mean_rank),]
  
  # full conparison table
  kable(social_value_sets, row.names = F) %>%
    kable_styling(full_width = F) %>% 
    add_header_above(c("State" = 1, "CAM" = 2, "unpac" = 2, "Difference" = 2))
  
  
# FIGURE 2: standard and UNPAc social value set comparison ----
  plot_df = data.frame(
    "healthstate" = all_states_str,
    "Standard" = val_set_mean,
    "standard_rank" = val_set_mean_rank,
    "UNPAc" = val_set_norm_mean
  ) 
  
  plot_df <- melt(plot_df, id.vars = c("standard_rank","healthstate"))
  
  plot_df <- plot_df[order(c(plot_df$standard_rank)),]
  plot_df$method <- plot_df$variable
  

  comp_labels = all_states_str[order(rank(-val_set_mean))][c(rep(c(T,rep(F,27)),8),T,rep(F,17),T)]
  comp_breaks =  c(1:243)[c(rep(c(T,rep(F,27)),8),T,rep(F,17),T)]
  
    plot_1 <- ggplot(plot_df) +
      geom_line(aes( x = standard_rank, y = value, col = method)) +
      theme_minimal() +
      ylab("Social health state value") +
      scale_x_continuous(
        name = "EQ-5D-3L health states, ordered by rank in standard model",
        breaks = comp_breaks, labels = comp_labels
        ) +
      theme(axis.text.x = element_text(angle = 45))
    
    plot_1
    ggsave(plot_1, filename = "./output/plot_1.png", bg = "white", width = 8, height = 6)
    
    
    
## INFLUENCE --------------------------

# ID's row indices for loop
  ID_rows = dat %>%
    group_by(ID) %>%
    dplyr::group_rows()
  names(ID_rows) = unique(dat$ID)
  util_range_by_ID = dat_norm$util_range[!duplicated(dat_norm$ID)]

# empty matrix to store relative and absolute social CAM and unpac values 
infl_mat = matrix(
  data=NA,
  ncol= 4,
  nrow = length(ID_rows),
  dimnames = list(
    names(ID_rows),
    c("mean_abs","mean_rel",
      "unpac_abs","unpac_rel")))

## *COMPUTING INDIVIDUALS INFLUENCE - LOOP 
  # this loop iteratively 
  # 1. removes an individuals from the data set
  # 2. re-fits the main  effects model 
  # 3. estimates the relative and absolute social value sets
  #    using the standard and UNPAc approach
  # 4. stores the results in the matrix

  for(i in 1:length(ID_rows)){
    cat("\r    ",  round(i/length(ID_rows),3)*100,"%",sep="") # displays status in loop
    
    x = ID_rows[[i]] # row index of the i'th participant 
    
    # 1 The standard approach
    fit_mean_i = lm(1-tto ~ MO2 + MO3 + SC2 + SC3 + UA2 + UA3 + PD2 + PD3 + AD2 + AD3 + N3,dat[-x,])
    val_set_mean_i = 1 - predict(fit_mean_i, newdata = all_states_df)
    val_set_mean_i[1] = 1
    val_set_rel_mean_i_diff = sum(abs(normaliseVS(val_set_mean_i) - val_set_rel_mean))
    val_set_mean_i_diff = sum(abs(val_set_mean_i - val_set_mean))
    
    # 2 the UNPAc approach
    fit_norm_i = lm(1-rel_val ~ MO2 + MO3 + SC2 + SC3 + UA2 + UA3 + PD2 + PD3 + AD2 + AD3 + N3, dat_norm[-x,])
    
    mean_anchoring_i = mean(dat_norm$util_range[-x])
    
    fit_norm_mean_i     = cbind(fit_norm_i$coefficients) * mean_anchoring_i
    val_set_norm_mean_i = predict.anchoredModel(fit_norm_mean_i, all_states_df)
    val_set_norm_mean_i[1] = 1
    val_set_rel_norm_mean_i_diff = sum(abs(normaliseVS(val_set_norm_mean_i) - val_set_rel_norm_mean))
    val_set_norm_mean_i_diff = sum(abs(val_set_norm_mean_i - val_set_norm_mean))
    
    infl_mat[i,] <- c(
      val_set_mean_i_diff, val_set_rel_mean_i_diff,
      val_set_norm_mean_i_diff, val_set_rel_norm_mean_i_diff
      )
    
  }
  
  write.csv(infl_mat, file = "./output/infl_mat.csv", row.names = F)
  

# FIGURE 3: Assess relationship between agents' utility ranges and their influences 
  infl_df <- data.frame(infl_mat)
  infl_df$util_range <- util_range_by_ID
  infl_df <- reshape2::melt(infl_df, id.vars = "util_range")
  infl_df$type = ifelse(grepl("rel",infl_df$variable), "rel", "abs")
  infl_df$method <- ifelse(infl_df$variable=="mean_abs" | infl_df$variable=="mean_rel","Standard","UNPAc")
  
  p1_e <- ggplot(infl_df[infl_df$type == "abs",]) +
    geom_point(aes(x = util_range, y = value, col = method), size = 0.1, alpha = 0.3) +
    geom_smooth(aes(x = util_range, y = value, col = method), se = F) +
    facet_wrap(~method) +
    theme_minimal() +
    xlab("Utility range") +
    ylab("Influence") +
    ggtitle("Influence on the absolute social value set") +
    theme(legend.position = "top") +
    # xlim(c(1,2)) +
    coord_cartesian(ylim= c(0,.1))
  
  
  p2_e<- ggplot(infl_df[infl_df$type == "rel",]) +
    geom_point(aes(x = util_range, y = value, col = method), size = 0.1, alpha = 0.3) +
    geom_smooth(aes(x = util_range, y = value, col = method), se = F) +
    facet_wrap(~method) +
    theme_minimal() +
    xlab("Utility range") +
    ylab("Influence") +
    ggtitle("Influence on the relative social value set") +
    theme(legend.position = "none") +
    # xlim(c(1,2)) +
    coord_cartesian(ylim= c(0,.075))
  
  pl_e <- get_legend(p1_e + theme(legend.position = "top"))
  
  
  plg_emp <- plot_grid(
    ncol = 1,
    rel_heights = c(0.45, 0.45,0.1),
    p1_e + theme(legend.position = "none"),
    p2_e + theme(legend.position = "none"),
    pl_e
    
   )


  save_plot(plg_emp,filename = "./output/plg_e.png",base_height = 5)



# # # Additional analyses for paper
#   #
#   # Assess differences in influences -----------
#   ggplot(infl_df, aes(value, fill = method)) +
#     geom_histogram( alpha = 0.5, bins = 100, position = "identity")  +
#     facet_wrap(~type)
#   # 
#   # # Correlation between influnces on social values in standard and unpac
#   # # relative social values
#   cor_rel_pear = cor(infl_mat[,1],infl_mat[,3], method = "pearson")
#   cor_rel_spear = cor(infl_mat[,1],infl_mat[,3], method = "spearman")
#   # # absolute social values
#   cor_abs_pear = cor(infl_mat[,2],infl_mat[,4], method = "pearson")
#   cor_abs_spear = cor(infl_mat[,2],infl_mat[,4], method = "spearman")
  