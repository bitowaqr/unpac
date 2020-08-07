
##----------------------------------------------------------
## Multi-step utility aggregation procedure
## Application using EQ-5D 3L data from the 1993 MVH Study 
##  
## Paul Schneider
##  University of Sheffield
##  p.schneider@sheffield.ac.uk
##  31.01.2020
##----------------------------------------------------------


## SETUP -----------------------------------------
#  rm(list=ls())   # clear global enviroment

# # decimal options
  options(scipen = 99)

# *load neccessary packages -----
  # library(tibble)
  # library(eq5d)
  library(reshape2)
  library(ggplot2)
  library(cowplot)
  # library(grid)
  # library(gridExtra)
  # library(corrplot)
  # library(knitr)
  # library(kableExtra)
  # library(ggExtra)
  # library(memisc)
  # library(Hmisc)
  library(RcppArmadillo)
  library(dplyr)


# *source convenient functions ----
  source("./R/support_f.R")

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

# . -------  
######## ZERO-DEVISION EXCLUSION   ----------
  # to avoid division by zero, participants with a zero utility range have to be excluded
  non_zero_range = aggregate( tto ~ ID, dat, function(x){ (max(x) - min(x)) > 0 } )
  cat("\n ***REMOVE",length(non_zero_range$ID[!non_zero_range$tto]),"IDs BECAUSE OF '0' RANGE!***")
  dat = filter(dat,ID %in% non_zero_range$ID[non_zero_range$tto])
#_______________________________________________

### CREATE AN EMPTY SOCIAL VALUE SET: A DATA FRAME WITH ALL 243 HEALTH STATES
  all_states_str = expand.grid(MO=1:3, SC=1:3,UA=1:3, PD=1:3, AD=1:3)
  all_states_str = apply(all_states_str, 1, function(x){paste(x,collapse="")})
  # head(all_states_str)  # health states as strings
  # new df for predictions
  all_states_df = hs_to_coef(all_states_str) # make new df from it
  # head(all_states_df)  # health states in parameter form
#_______________________________________________


## MODEL FITTING ------------------------------- 
# *FIT CAM -----
  # Dolan's main effects model
  fit_CAM = lm(1-tto ~ MO2 + MO3 + SC2 + SC3 + UA2 + UA3 + PD2 + PD3 + AD2 + AD3 + N3,dat)
  # summary(fit_CAM) 
  
  # Estimate CAM social value set
  val_set_CAM = 1 - predict(fit_CAM,newdata = all_states_df) 
  val_set_CAM[1] = 1 # full health = 1 correction

### *FIT MUAP --------
  # **Step 1: Individual min-max utility normalisation ----
  dat = dat %>%
    group_by(ID) %>%
    mutate(
      util_range = max(tto)-min(tto),
      rel_val = (tto - min(tto)) / (max(tto) - min(tto)) 
    )

  # ** Step 2: Relative social tariff ---------
  fit_MUAP = lm(1-rel_val ~ MO2 + MO3 + SC2 + SC3 + UA2 + UA3 + PD2 + PD3 + AD2 + AD3 + N3,dat)
  # summary(fit_MUAP)
  # Estimate relative social values for all states
  rel_val_set_MUAP = 1-predict(fit_MUAP,newdata = all_states_df)
  
  # ** Step 3: Social scaling factor -------------
  scaling_f = mean(dat$util_range)

  # ** Step 4: Rescaling -------------
  fitted_MUAP = rel_val_set_MUAP * scaling_f  + (1-scaling_f)
  fitted_MUAP[1] = 1 # full health correction

#####---- Goodness of fit R2 stats  
  fitted_MUAP_hat = (1-predict(fit_MUAP)) * scaling_f  + (1-scaling_f)
  # cor(predict(fit_MUAP),dat$rel_val)^2 # R2 of rel MUAP social tariff
  # cor(fitted_MUAP_hat,dat$tto)^2       # R2 of abs MUAP social tariff
  # cor(predict(fit_CAM),dat$tto)^2      # R2 of abs CAM social tariff
  
  

## RESCALE MUAP MODEL COEFFICIENTS and CREATE COMPARISON TABLE --------
  MUAP_resc = cbind(fit_MUAP$coefficients,confint(fit_MUAP)[,1],confint(fit_MUAP)[,2]) * scaling_f
  MUAP_resc = formatC(MUAP_resc,digits=2,format="f")  
  MUAP_resc = paste(MUAP_resc[,1]," (",MUAP_resc[,2],"-",MUAP_resc[,3],")",sep="")
  MUAP_resc = c(MUAP_resc,"","")
  b_ratios = (fit_MUAP$coefficients*scaling_f)/fit_CAM$coefficients
  b_ratios = c(formatC(b_ratios, digits=2,format="f"),"","")
  fit_tbl = cbind(
    make_tbl(fit_CAM,get.rownames = T,give.colname = "CAM"),
    make_tbl(fit_MUAP,get.rownames = F,give.colname = "MUAP (normalised)"),
    "MUAP (rescaled)" = MUAP_resc,
    "Beta ratios" = b_ratios)
  # fit_tbl     
  write.csv(fit_tbl,file = "./output/model_comparison.csv",row.names = F)

  
## COMPARE CAM and MUAP SOCIAL VALUE SETS -----------
  social_value_sets = data.frame(healthstate = all_states_str,
                        CAM_value = val_set_CAM,
                        MUAP_value = fitted_MUAP,
                        value_diff = fitted_MUAP-val_set_CAM,
                        CAM_rank = rank(-val_set_CAM),
                        MUAP_rank = rank(-fitted_MUAP),
                        rank_diff = rank(-fitted_MUAP)-rank(-val_set_CAM)
  ) 
  social_value_sets = arrange(social_value_sets,CAM_rank)
  # social_value_sets
  write.csv(social_value_sets,file = "./output/social_value_sets.csv",row.names = F)


  
# *DESCRITPTIVE STATS: SOCIAL VALUE SETS ----------------
  descriptive_stats = make_sts(CAM= social_value_sets$CAM_value,MUAP= social_value_sets$MUAP_value,hs_ref = social_value_sets$healthstate)
  # t(descriptive_stats)
  
  
# *PLOT: VISUAL COMPARISON OF SOCIAL VALUE SETS --------
  comparison_plot_df = social_value_sets[,c("healthstate","CAM_value","MUAP_value","CAM_rank")]
  comp_labels = comparison_plot_df$healthstate[c(rep(c(T,rep(F,9)),24),F,F,T)]
  comp_breaks =  comparison_plot_df$CAM_rank[c(rep(c(T,rep(F,9)),24),F,F,T)]
  comparison_plot_df_long = melt(comparison_plot_df,id.vars = c("healthstate","CAM_rank"))
  
  comparison_plot = 
    ggplot(comparison_plot_df_long) +
    geom_hline(yintercept = 0,size=0.2,linetype="dashed")+
    geom_line(aes(x=CAM_rank,y=value,col=variable),size=0.5,alpha=1) +
    geom_point(aes(x=CAM_rank,y=value,col=variable),size=1,alpha=0.5,shape=20) +
    theme_minimal()+
    ylab("Social value") +
    scale_x_continuous(breaks = comp_breaks,minor_breaks = NULL,
                       name="EQ-5D 3L health state (ordered highest to lowest value)",
                       labels = comp_labels) +
    theme(axis.text.x = element_text(angle = 35,size = 8),
          legend.position = "top",
          legend.text = element_text(size = 12)) +
    scale_color_manual(values = c("#396AB1","#E1974C"),name="Utility aggregation procedure",labels=c("CAM","MUAP")) +
    NULL
  
  ggsave(filename = "./output/comparisonplot.png",comparison_plot,
         height =5.25, width=7)


# . ---------------------------------
## INFLUENCE --------------------------

# *Influence Setup ---------
  # compute IDs influence on model
  # re-fit CAM model
  fit_CAM = lm(formula = 1 - tto ~ MO2 + MO3 + SC2 + SC3 + UA2 + UA3 + PD2 + PD3 + AD2 + AD3 + N3, data = dat)
  CAM_tarif = 1-predict(fit_CAM,all_states_df)  
  CAM_tarif[1] = 1 # set full health = 1
  CAM_tarif_rel = (CAM_tarif - min(CAM_tarif)) / (max(CAM_tarif)-min(CAM_tarif)) # relative values
  
  # re-fit MUAP model
  fit_MUAP = lm(formula = 1-rel_val ~ MO2 + MO3 + SC2 + SC3 + UA2 + UA3 + PD2 + PD3 + AD2 + AD3 + N3, data = dat)
  MUAP_tarif_rel = 1-predict(fit_MUAP,all_states_df)
  MUAP_tarif_rel[1] = 1
  scaling_f = mean(dat$util_range)
  MUAP_tarif = (MUAP_tarif_rel * scaling_f  + (1-scaling_f))
  
  # ID's row indices for loop
  ID_rows = dat %>%
    group_by(ID) %>%
    dplyr::group_rows()
  names(ID_rows) = unique(dat$ID)
  
  util_range_by_ID = dat$util_range[!duplicated(dat$ID)]
  
  fast_df = as.matrix(cbind(1,dat[,4:14]))     # model matrix to efficiently refit model
  fast_hs = as.matrix(cbind(1,all_states_df))  # prediction matrix to efficiently predict social values
  
  # empty matrix to store relative and absolute social CAM and MUAP values 
  CAM_abs_M = CAM_rel_M = MUAP_abs_M = MUAP_rel_M = matrix(data=NA,
                                                         ncol=length(unique(dat$ID)),
                                                         nrow = nrow(all_states_df),
                                                         dimnames = list(c(all_states_str),unique(dat$ID)))
  
## *COMPUTING INDIVIDUALS INFLUENCE - LOOP -------
  # this loop iteratively 
  # 1. removes an individuals from the data set
  # 2. re-fits the main  effects model 
  # 3. estimates the relative and absolute social value sets
  #    using the CAM and MUAP approach
  # 4. stores the results in the matrix
  
  for(i in 1:length(ID_rows)){
    cat("\r    ",  round(i/length(ID_rows),3)*100,"%",sep="") # displays status in loop
    x = ID_rows[[i]] # row index of the i'th participant 
    
    # CAM
    # fast main effecs model fitting
    temp.lm.coef = fastLm(y= 1-dat$tto[-x],X= fast_df[-x,])$coefficients 
    # abs social value set estimation
    CAM_abs_M[,i] = 1- (fast_hs %*% temp.lm.coef)
    CAM_abs_M[1,i] = 1 
    # relative value set estimation
    CAM_rel_M[,i] = (CAM_abs_M[,i] - min(CAM_abs_M[,i])) / (max(CAM_abs_M[,i])-min(CAM_abs_M[,i]))
    
    # MUAP
    # fast main effecs model fitting
    temp.lm.coef2 = fastLm(y= 1-dat$rel_val[-x],X= fast_df[-x,])$coefficients
    # relative value set estimation
    MUAP_rel_M[,i] = 1- (fast_hs %*% temp.lm.coef2)
    MUAP_rel_M[1,i] = 1
    # retrieve scaling factor without individual i in the dataset
    temp.scaling_f = mean(dat$util_range[-x])
    # abs social value set estimation
    MUAP_abs_M[,i] = (MUAP_rel_M[,i]) * temp.scaling_f  + (1-temp.scaling_f)
  }
  


# *Analyze influences -----------

  # compute individual influence metric as
  # the sum of absolute differences between social values with and without individual i
  influence_df = data.frame(
    ID = names(ID_rows),
    CAM_abs_diff = get_diff(CAM_abs_M,CAM_tarif),
    CAM_rel_diff = get_diff(CAM_rel_M,CAM_tarif_rel),
    MUAP_abs_diff = get_diff(MUAP_abs_M,MUAP_tarif),
    MUAP_rel_diff = get_diff(MUAP_rel_M,MUAP_tarif_rel),
    scaling_f = util_range_by_ID
  )

  
# # *Show an example: what does 'influence' capture? -------
  # sample_plot_df = data.frame(hs = all_states_str,
  #                             CAM_tarif = CAM_tarif,
  #                             order = rank(-CAM_tarif),
  #                             CAM_abs_M[,"1473"])
  # sample_plot_df = melt(sample_plot_df,id.vars = c("hs","order"))
  # ggplot(sample_plot_df) +
  #   geom_point(aes(x=order,y=value,col=variable),size=0.5)+
  #   geom_line(aes(x=order,y=value,col=variable),size=0.5) +
  #   scale_color_manual(values = 2:3,labels=c("With individual","Without individual"),name="Social tariff")+
  #   scale_x_continuous(name="Health state",breaks = sample_plot_df$order,labels = sample_plot_df$hs,limits = c(50,55)) +
  #   ylim(0.32,0.351)  +
  #   theme(axis.text.x = element_text(angle=45,hjust = 1))

# *Utility ranges descriptive stats ------
  range_stats = data.frame(
    Mean = formatC(mean(influence_df$scaling_f),digits=2,format="f"),
    SD = formatC(sd(influence_df$scaling_f),digits=2,format="f"),
    Median = formatC(median(influence_df$scaling_f),digits=2,format="f"),
    IQR = formatC(IQR(influence_df$scaling_f),digits=2,format="f"),
    MIN = formatC(min(influence_df$scaling_f),digits=2,format="f"),
    MAX = formatC(max(influence_df$scaling_f),digits=2,format="f")
  )
  # range_stats

  
# dispersion descriptive stats ----------
  disp_stats = get_dispersion(MUAP = influence_df$MUAP_rel_diff,
                              CAM = influence_df$CAM_rel_diff)
  # disp_stats
  write.csv( disp_stats,file="./output/disp_stats.csv")


# Correlation between influnces on relative social values in CAM and MUAP
  cor_cm_muap_pearson = formatC(cor(influence_df$CAM_rel_diff,influence_df$MUAP_rel_diff,method = "pearson"),digits=2,format="f")
  cor_cm_muap_spearma = formatC(cor(influence_df$CAM_rel_diff,influence_df$MUAP_rel_diff,method = "spearman"),digits=2,format="f")


#### ASSOCIATION BETWEEN UTILITY RANGE AND INFLUENCE -----------------
  # *fitting polynomial models ------
  influence_plot_df = melt(influence_df,id.vars=c("ID","scaling_f"))
  r_fit_cam_1 = lm(CAM_rel_diff~ poly(scaling_f,1),influence_df )
  r_fit_cam_2 = lm(CAM_rel_diff~ poly(scaling_f,2),influence_df )
  r_fit_cam_3 = lm(CAM_rel_diff~ poly(scaling_f,3),influence_df )
  r_fit_muap_1 = lm(MUAP_rel_diff ~ poly(scaling_f,1),influence_df )
  r_fit_muap_2 = lm(MUAP_rel_diff ~ poly(scaling_f,2),influence_df )
  r_fit_muap_3 = lm(MUAP_rel_diff ~ poly(scaling_f,3),influence_df )
  
  mtbl = Reduce(function(x, y) cbind(x, y), list(make_mtbl(r_fit_cam_1), 
                                                 make_mtbl(r_fit_cam_2), 
                                                 make_mtbl(r_fit_cam_3),
                                                 make_mtbl(r_fit_muap_1), 
                                                 make_mtbl(r_fit_muap_2), 
                                                 make_mtbl(r_fit_muap_3)))
  rownames(mtbl)[1:4] = c("Intercept",
                          "Range",
                          "Range$^2$",
                          "Range$^3$")
  names(mtbl) = c("CM linear","CM square","CM cubic","MUAP linear","MUAP sqaure","MUAP cubic")
  
  #  mtbl
  write.csv(mtbl,file = "./output/poly_influence_models.csv")

  
  # *Visual comparison ------------
  p_CAM_df = influence_plot_df[grepl("CAM_rel",influence_plot_df$variable),]
  
  p_CAM =
    ggplot(p_CAM_df)+
    geom_point(aes(x=(scaling_f),y=value,col="1"),alpha=0.5,size=0.5) +
    # geom_hline(color="black",yintercept = mean(influence_plot_df[grepl("CM_rel",influence_plot_df$variable),]$value),linetype=2) +
    geom_smooth(aes(x=(scaling_f),y=value,col="2"),se=F,method = "lm",linetype=1,size=0.9)+
    geom_smooth(aes(x=(scaling_f),y=value,col="3"),se=F,method = "lm",formula = "y~poly(x,2)",size=0.9)+
    geom_smooth(aes(x=(scaling_f),y=value,col="4"),se=F,method = "lm",formula = "y~poly(x,3)",size=0.9)+
    
    annotate(geom="text",x=0.2,y=0.067,label=c("CAM") ) +
    coord_cartesian(ylim=c(0,0.07)) +
    ylab("Influence on relative social values") + 
    xlab("Range of utility values") +
    scale_color_manual(values=c(RColorBrewer::brewer.pal(4,"Dark2")),
                       labels = c("","Linear model","Quadratic model","Cubic model"),
                       name="") +
    theme_minimal() +
    theme(legend.position = "none")  +
    NULL
  
  p_muap_df = influence_plot_df[grepl("MUAP_rel",influence_plot_df$variable),]
  p_MUAP =
    ggplot(p_muap_df)+
    geom_point(aes(x=(scaling_f),y=value,col="1"),alpha=0.5,size=0.5) +
    # geom_hline(color="black",yintercept = mean(influence_plot_df[grepl("MUAP_rel",influence_plot_df$variable),]$value),linetype=2) +
    geom_smooth(aes(x=(scaling_f),y=value,col="2"),se=F,method = "lm",linetype=1,size=0.9)+
    geom_smooth(aes(x=(scaling_f),y=value,col="3"),se=F,method = "lm",formula = "y~poly(x,2)",size=0.9)+
    geom_smooth(aes(x=(scaling_f),y=value,col="4"),se=F,method = "lm",formula = "y~poly(x,3)",size=0.9)+
    
    annotate(geom="text",x=0.2,y=0.067,label=c("MUAP")) +
    coord_cartesian(ylim=c(0,0.07)) +
    ylab("Influence on relative social values") + 
    xlab("Range of utility values") +
    scale_color_manual(values=c(RColorBrewer::brewer.pal(4,"Dark2")),
                       labels = c("","Linear model","Quadratic model","Cubic model"),
                       name="") +
    theme_minimal() +
    theme(legend.position = "none",axis.title.y = element_text(color="white"))  +
    NULL
  
  plot_legend_top =   get_legend( p_MUAP + 
                                    theme(legend.position = "top",
                                          legend.text = element_text(size=12)
                                          ) +
                                    guides(line = guide_legend(
                                      override.aes = list(size = 5)
                                      )))
  
  multiplot = plot_grid(
    plot_grid(p_CAM,p_MUAP,nrow=1,align = "v"),
    plot_legend_top,nrow = 2,rel_heights = c(9,1))
  
  # multiplot
  ggsave(filename = "./output/multiplot.png",plot=multiplot,height = 4,width=12)

##### Additional analyses
  # number of individuals with <1 utility range
  uniq_id_len = length(unique(dat$ID))
  n_r_1 = sum(influence_df$r<1)
  perc_r_1 = round(sum(influence_df$r<1)/length(influence_df$r),2)*100

  # # irrationals: '33333' not lowest value ---------
  irrationals = group_by(dat,ID) %>%
    summarise(irrational = "33333" %in% healthstate[tto != min(tto)]) %>%
    pull()
  # sum(irrationals);sum(irrationals)/length(irrationals)
  
  
  ## Additional figure for paper
  ## Illustration of utility range differences between participants
  sample_profs = dat %>% 
    filter(healthstate %in% c("11111","21111","22222","33333")) %>%
    group_by(ID) %>%
    filter(length(ID) ==4) %>%
    select(1:3) %>% 
    mutate(hu = tto[healthstate == "33333"]) %>%
    filter(ID %in% c(5075,377,3895,5651,439,729,4212,129))
  
  
  p_sample_profs = 
    ggplot(sample_profs) +
    geom_hline(yintercept = 0,size=0.2,linetype="dashed") +
    geom_point(aes(x=healthstate,y=tto,group=as.factor(ID),col=hu))+
    geom_line(aes(x=healthstate,y=tto,group=as.factor(ID),col=hu)) +
    scale_x_discrete(labels = c("11111\n(full health)","21111\n(mild)","22222\n(moderate)","33333\n(severe)")) +
    theme_minimal() +
    xlab("EQ-5D 3L health state") +
    ylab("TTO utility value") +
    theme(axis.text.x = element_text(angle = 35,size = 8),
          legend.position = "none",
          legend.text = element_text(size = 12))  +
    NULL
  ggsave(filename = "./output/p_sample_profs.png",p_sample_profs,height=4,width = 8)  
  
  
### 'A simple example' (ase)
  # plot and figure for paper
  Alice = c(1,0.966666666,0.9333333,0.9)
  Bob = c(1,-0.333333333,0.333333333,-1)
  ase_df = (t(cbind(Alice,Bob)))
  colnames(ase_df) = c("h_full","h_a","h_b","h_O")
  # CAM social value set
  ase_CAM = colMeans(ase_df)
  # normalised utilities
  ase_norm_utils = data.frame(t(apply(ase_df,1,function(x){(x-min(x))/(max(x)-min(x))})))
  # relative social value set
  ase_rel_MUAP = colMeans(ase_norm_utils)
  # utility ranges
  ase_util_ranges = unlist(apply(ase_df,1,function(x){max(x)-min(x)}))
  # scaling factor
  ase_scaling_f = mean(ase_util_ranges)
  # MUAP social value set
  ase_MUAP = ase_rel_MUAP * ase_scaling_f + (1-ase_scaling_f)
  
  p_ase = 
    ggplot() +
    geom_hline(yintercept = 0,size=0.3) +
    geom_point(aes(x=1:4,y=ase_df[1,],col="Alice"),size=2) +
    geom_path(aes(x=1:4,y=ase_df[1,],col="Alice",linetype="Alice"),size=1) +
    geom_point(aes(x=1:4,y=ase_df[2,],col="Bob"),size=2) +
    geom_path(aes(x=1:4,y=ase_df[2,],col="Bob",linetype="Bob"),size=1) +
    geom_point(aes(x=1:4,y=ase_CAM,col="CAM"),size=2) +
    geom_path(aes(x=1:4,y=ase_CAM,linetype="CAM",col="CAM"),size=1.1) +
    geom_point(aes(x=1:4,y=ase_MUAP,col="MUAP"),size=2) +
    geom_path(aes(x=1:4,y=ase_MUAP,linetype="MUAP",col="MUAP"),size=1.1) +
    scale_x_continuous(labels = colnames(ase_df), name="Health state") +
    ylab("Utility value")  +
    theme_minimal() + 
    scale_linetype_discrete(guide = FALSE) +
    scale_color_brewer(palette = "Dark2",name = "") 
  
  ggsave(filename = "./output/p_ase.png",plot=p_ase,height = 5,width=7)
  
  
  
## fin.
