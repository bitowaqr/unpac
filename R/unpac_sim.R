
# -----------------------------------
# UNPAc utility aggregation procedure
# by Paul Schneider
# University of Sheffield
# MIT license
# July 2021
# -----------------------------------


# SETUP  ---------

# # load required packages 
  library(ggplot2)
  library(cowplot)
  library(truncnorm)
  library(kableExtra)
  library(dplyr)
  
  options(scipen = 99)

# # clear global enviroment space 
  # rm(list=ls())

# # random seed 
  set.seed(2021)

# # utility function to normalise preferences
  normaliseVS <- function(mat){
    1 -  ((1- mat) / (max(mat)- min(mat)))
  }


# SIMULATION PARAMETERS --------

  hs <- c("FH","P","M","PM") # health states: FH, M ,P, PM
  n_hs = length(hs)         # number of states
  n_agents <- 10000         # number of agents to simulate
  prop_Gm <- 0.45           # proportion of agents in group Gm
  n_Gm <- round(prop_Gm * n_agents) # number of agents in group Gm
  n_Gp <- n_agents - n_Gm           # agents in group Gp


# CREATE PREFERENCE FUNCTIONS -------------------------

# # Gm: FH > M > P > PM
  Gm_pref <- matrix(NA, nrow = n_Gm, ncol = n_hs)
  colnames(Gm_pref) <- hs
  Gm_pref[,1] <- 1 # u(FH) = 1
  Gm_pref[,n_hs] <- rtruncnorm(n = n_Gm, a = -1, b = 1, mean = -0.9, sd = 0.3) # u(PM)
  Gm_pref_mid <- cbind(replicate(n = n_hs - 2,runif(n_Gm, Gm_pref[,n_hs], 1))) # u(M) and u(P)
  Gm_pref_mid <- t(apply(Gm_pref_mid,1,function(x)x[order(x, decreasing = F)]))
  Gm_pref[,-c(1,n_hs)] <- Gm_pref_mid

# # Gp: FH > P > M > PM
  Gp_pref <- matrix(NA, nrow = n_Gp, ncol = n_hs)
  colnames(Gp_pref) <- hs
  Gp_pref[,1] <- 1 # u(FH) = 1
  Gp_pref[,n_hs] <- rtruncnorm(n = n_Gp, a = -1, b = 1, mean = -0.25, sd = 0.3) # u(PM)
  Gp_pref_mid <- cbind(replicate(n = n_hs - 2, runif(n_Gp, Gp_pref[,n_hs], 1))) # u(M) and u(P)
  Gp_pref_mid <- t(apply(Gp_pref_mid,1,function(x)x[order(x, decreasing = T)]))
  Gp_pref[,-c(1,n_hs)] <- Gp_pref_mid

# # merge Gm and Gp (= society)
  pref_mat <- rbind(Gm_pref, Gp_pref)

  
  
  
  
# STANDARD UTILITY AGGREGATION APPROACH --------------

  standard_vset     <- apply(pref_mat, 2, mean)
  standard_relative_vset <-  normaliseVS(standard_vset)
  standard_hs_ranks <- rank(-standard_vset)

  
# UNPAc APPROACH -----------------

  # Step 1: normalise
  normalised_preferences <- t(apply(pref_mat,1,normaliseVS))
  # Step 2: aggregate (mean)
  normalised_vset <- apply(normalised_preferences, 2, mean)
  # Step 3a: retrieve individual anchor points 
  anchor <- apply(pref_mat,1,min)
  # step 3b: aggregate anchor points and re-anchor normalised social value set
  unpac_vset <- 1 + normalised_vset * (1 - mean(anchor)) - (1- mean(anchor))
  unpac_relative_vset <-  normalised_vset
  unpac_hs_ranks <- rank(-unpac_vset)


  
    
# ASSESS INFLUENCE -------
  # agents' influences on absolute and relative social value set --------

  infl_standard <- infl_standard_rel <- 
    infl_norm_mean <- infl_norm_mean_rel <- rep(NA, n_agents)
  
  # N.B.: This loop may take quite a while(!) to run for 10,000 agents
  for(i in 1:n_agents){
    
    cat("\r  Status: ", i," of", n_agents, "     ")
    
    # standard
    standard_vset_i <- apply(pref_mat[-i,], 2, mean)
    
    infl_standard[i]  <- sum(abs(standard_vset_i - standard_vset))
    infl_standard_rel[i] <-  sum(abs(normaliseVS(standard_vset_i) - standard_relative_vset))
    
    # UNPAc
    normalised_preferences_i <-t(apply(pref_mat[-i,],1,normaliseVS)) 
    normalised_vset_i <- apply(normalised_preferences_i, 2, mean)
    anchor_i <- mean(anchor[-i])
    unpac_vset_i <- 1 + normalised_vset_i * (1 - anchor_i) - (1- anchor_i)
    
    infl_norm_mean[i]  <- sum(abs(unpac_vset_i - unpac_vset))
    infl_norm_mean_rel[i] <- sum(abs(normalised_vset_i - unpac_relative_vset))
    
  }



#############################
### OUTPUTS ------------
##################


# TABLE 1 ---------
  tbl1a <- apply(Gm_pref,2,function(x) formatC(c(mean = mean(x), sd = sd(x)),digits = 2, format = "f"))
  tbl1a_rel <- formatC(normaliseVS(apply(Gm_pref,2,mean)),digits = 2, format = "f")
  tbl1b <- apply(Gp_pref,2,function(x) formatC(c(mean = mean(x), sd = sd(x)),digits = 2, format = "f"))
  tbl1b_rel <- formatC(normaliseVS(apply(Gp_pref,2,mean)),digits = 2, format = "f")
  tbl1c <- apply(pref_mat,2,function(x) formatC(c(mean = mean(x), sd = sd(x)),digits = 2, format = "f"))
  tbl1c_rel <- formatC(normaliseVS(apply(pref_mat,2,mean)),digits = 2, format = "f")
  tbl1 <- cbind(
    cbind(paste0(tbl1a[1,], " (",tbl1a[2,],")"), tbl1a_rel),
    cbind(paste0(tbl1b[1,], " (",tbl1b[2,],")"),  tbl1b_rel),
    cbind(paste0(tbl1c[1,], " (",tbl1c[2,],")"),  tbl1c_rel)
  )
  rownames(tbl1) <- c("fh","i","j","pit")
  colnames(tbl1) <- c("state","GM mean","GM norm","GP mean","GP norm","All norm")
  
  kable(tbl1, col.names = rep(c("Mean (SD)","Normalised value"),3)) %>%
    kable_styling() %>%
    add_header_above(header = c("state","Gm (n = 4,500)" = 2, "Gp (n = 5,500)" = 2, "all" = 2))


  
# TABLE 2 ---------
  
  tbl2 <- cbind(
    "rank" = standard_hs_ranks,
    "Absolute (relative) value" = paste0(
      formatC(standard_vset,digits = 2, format = "f"),
      " (", 
      formatC(standard_relative_vset,digits = 2, format = "f"),
      ")"),
    "rank" = unpac_hs_ranks,
    "Absolute (relative) value" = paste0(
      formatC(unpac_vset,digits = 2, format = "f"),
      " (", 
      formatC(normalised_vset,digits = 2, format = "f"),
      ")"),
    "rank" = formatC(unpac_hs_ranks - standard_hs_ranks, digits = 0, format = "f", flag = "+"),
    "Absolute (relative) value" = paste0(
      formatC(unpac_vset - standard_vset, digits = 2, format = "f", flag = "+"),
      " (", formatC(normalised_vset - standard_relative_vset, digits = 0, format = "f", flag = "+"),
      ")"
    )
  )
  
  rownames(tbl2) <- hs
  kable(tbl2) %>%
    kable_styling(full_width = F) %>% 
    add_header_above(c("State" = 1, "Standard" = 2, "UNPAc" = 2, "Difference" = 2))
  
  
# Relationship between agents utility ranges and their influences ----------------
  
  infl_df <- data.frame(
    agents_ranges = 1-apply(pref_mat, 1, min),
    value = c(
      infl_standard,infl_standard_rel,
      infl_norm_mean, infl_norm_mean_rel),
    method = rep(c("standard",
                   "UNPAc"), each = n_agents*2),
    type = rep(c("abs", "rel"), each = n_agents)
  )
  
  
# simple summary stats --------
  
  # summary stats utility function
  sumaryStats <- function(x){
    cbind(
      "mean" = mean(x),
      "sd" = sd(x),
      t(quantile(x, c(0.025, 0.5, 0.975)))
    )
  }
  
  # summary stats
  sumaryStats(infl_df$value[infl_df$type == "abs" & infl_df$method == "standard"])
  sumaryStats(infl_df$value[infl_df$type == "abs" & infl_df$method == "UNPAc"])
  
  sumaryStats(infl_df$value[infl_df$type == "rel" & infl_df$method == "standard"])
  sumaryStats(infl_df$value[infl_df$type == "rel" & infl_df$method == "UNPAc"])
  
  
  
# # Figure 1: influence on absolute and relative value sets in the standard and unpac approach  ------
  
  p1 <- ggplot(infl_df[infl_df$type == "abs",]) +
    geom_point(aes(x = agents_ranges, y = value, col = method), size = 0.1, alpha = 0.3) +
    geom_smooth(aes(x = agents_ranges, y = value, col = method), se = F) +
    facet_wrap(~method) +
    theme_minimal() +
    xlab("Utility range") +
    ylab("Influence") +
    ggtitle("Influence on the absolute social value set")
  
  p2 <- ggplot(infl_df[infl_df$type == "rel",]) +
    geom_point(aes(x = agents_ranges, y = value, col = method), size = 0.1, alpha = 0.3) +
    geom_smooth(aes(x = agents_ranges, y = value, col = method), se = F) +
    facet_wrap(~method) +
    theme_minimal() +
    xlab("Utility range") +
    ylab("Influence") +
    ggtitle("Influence on the relative social value set")
  
  pl <- get_legend(p1 + theme(legend.position = "top"))
  
  plg <- plot_grid(
    ncol = 1,
    rel_heights = c(0.45, 0.45,0.1),
    p1 + theme(legend.position = "none"),
    p2 + theme(legend.position = "none"),
    pl
  )
  
  save_plot(plg,filename = "./output/plg.png",base_height = 5)
  
  
  
  
  
  
###########################################################
# Additional Figures (not reported in the paper) -----
###############################################
  
# # u(PM) histogramm -------------
  # # u(PM) distribution
  # ggplot() +
  #   geom_histogram(aes(as.numeric(Gm_pref[,n_hs]), col = "Gm", fill = "Gm"), alpha = 0.5) +
  #   geom_histogram(aes(as.numeric(Gp_pref[,n_hs]), col = "Gp", fill = "Gp"), alpha = 0.5)  +
  #   theme_minimal()
  
# #within group plots -----------    
  # within_group_vs <- data.frame(
  #   value = c(apply(Gm_pref,2,mean),
  #             apply(Gp_pref,2,mean),
  #             apply(pref_mat,2,mean)
  #   ),
  #   group = rep(c("Gm", "Gp","Social"), each = n_hs),
  #   hs_x = 1:n_hs
  # )
  # 
  # ggplot(within_group_vs) +
  #   geom_point(aes(x = hs_x, y= value, col = group)) +
  #   geom_line(aes(x = hs_x, y= value, col = group), alpha = 1) +
  #   scale_x_continuous(breaks = 1:n_hs, labels = hs) +
  #   theme_minimal() +
  #   theme(legend.position = "top")

#  # within_group_rel_vs <- data.frame(
  #   value = c( normaliseVS(apply(Gm_pref,2,mean)),
  #              normaliseVS(apply(Gp_pref,2,mean)),
  #              normaliseVS(apply(pref_mat,2,mean))
  #   ),
  #   group = rep(c("Gm", "Gp","Social"), each = n_hs),
  #   hs_x = 1:n_hs
  # )
  # 
  
# # alternative comparison
  # vs_df <- data.frame(
  #   value = c(standard_vset, standard_relative_vset,
  #             # median_vs, median_rel_vs,
  #             # norm_median_vs, norm_median_rel_vs,
  #             unpac_vset, unpac_relative_vset),
  #   method = rep(c("Standard",
  #                  "UNPAc"),
  #                each = n_hs*2),
  #   type = rep(c("absolute value set","relative value set"), each = n_hs),
  #   hs = hs,
  #   hs_x = 1:n_hs
  # )
  # 
  # fiGp <- ggplot(vs_df[vs_df$type =="absolute value set",]) +
  #   geom_point(aes(x = hs_x, y= value, col = method)) +
  #   geom_line(aes(x = hs_x, y= value, col = method, linetype = method), alpha = 0.5) +
  #   # facet_wrap(~type, nrow = 1) +
  #   scale_x_continuous(breaks = 1:n_hs, labels = hs) +
  #   theme_minimal() +
  #   theme(legend.position = "top")
  # 
  # fiGp
