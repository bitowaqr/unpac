
# some convenient functions to extract results 

  # converts health state strings into paramter form 
  hs_to_coef = function(hs){
    df = data.frame(hs = hs, d1=NA,d2=NA,d3=NA,d4=NA,d5=NA,stringsAsFactors = F)
    
    df = df %>% 
      mutate(d1=substr(hs,1,1),
             d2 = substr(hs,2,2),
             d3 = substr(hs,3,3),
             d4 = substr(hs,4,4),
             d5 = substr(hs,5,5)) %>%
      mutate(MO2 = d1>1, MO3 = d1>2,
             SC2 = d2>1, SC3 = d2>2,
             UA2 = d3>1, UA3 = d3>2,
             PD2 = d4>1, PD3 = d4>2,
             AD2 = d5>1, AD3 = d5>2) %>%
      mutate(N3 = d1 == 3 | d2 == 3 | d3  == 3 | d4 == 3 | d5 == 3 )
    df = df[,-c(1:6)]
    df
  }
  
  
  # # utility function to normalise preferences
  normaliseVS <- function(mat){
    1 -  ((1- mat) / (max(mat)- min(mat)))
  }
  
  # predict values for anchored model
  predict.anchoredModel <- function(fit, hs){
    
    if(is.null(dim(hs))){
      return(1-sum(fit[,1] * c(1,as.numeric(hs))))
    }
    
    apply(hs,1,function(x){
      1-sum(fit[,1] * c(1,as.numeric(x)))
    })
  }
  
  # retrieves some useful stats from social value sets
  make_sts = function(CAM,MUAP,hs_ref){
    abs.diff = abs(MUAP-CAM)
    MUAP_rank = rank(-MUAP)
    CAM_rank = rank(-CAM)
    err=c()
    fitted_norm_order = MUAP[order(MUAP)]
    for(i in 2:length(fitted_norm_order)){
      err = c(err,abs(fitted_norm_order[i] - fitted_norm_order[i-1]))
    }
    err.m = mean(err)
    err.sd = sd(err)
    
    sts_df = data.frame(
      mean_abs_error = mean(abs.diff),
      max_abs_error = max(abs.diff),
      sd_abs_error = sd(abs.diff),
      mean_error = mean(MUAP-CAM),
      sd_error = sd(MUAP-CAM),
      comp_mean_abs_err = err.m,
      comp_sd_abs_err = err.sd,
      stringsAsFactors = F)
    sts_df = apply(sts_df,2,function(x){formatC(x,digits = 3,format = "f")})
    max_abs_error_hs = hs_ref[which(abs.diff == max(abs.diff))]
    CAM_at_max_diff = formatC(CAM[hs_ref==max_abs_error_hs],digits=3,format="f")
    MUAP_at_max_diff = formatC(MUAP[hs_ref==max_abs_error_hs],digits=3,format="f")
    
    rank_diff = MUAP_rank-CAM_rank
    rank_diff_n = sum((rank_diff) !=0)
    min_rank_diff = min(rank_diff)
    min_rank_diff_hs = all_states_str[which(rank_diff== min(rank_diff))][1]
    max_rank_diff = max(rank_diff)
    max_rank_diff_hs = all_states_str[which(rank_diff== max(rank_diff))][1]
    sts_df = data.frame(t(sts_df),
                        max_abs_error_hs = as.character(max_abs_error_hs),
                        CAM_at_max_diff = CAM_at_max_diff,
                        MUAP_at_max_diff =MUAP_at_max_diff,
                        rank_diff_n = rank_diff_n,
                        rank_diff_perc = paste(round(100*(rank_diff_n/length(abs.diff)),1)),
                        max_rank_diff = max_rank_diff,
                        max_rank_diff_hs = max_rank_diff_hs,
                        min_rank_diff = min_rank_diff,
                        min_rank_diff_hs = min_rank_diff_hs,
                        stringsAsFactors = F
    )
    sts_df
  }
  
  
  
  # # computes individuals' influences
  get_diff = function(M_i,      # matrix with social tariffs estimates without individual i
                      tariff,   # reference social value set
                      margin=2, # 2=columns represent tariffs without individul i
                      scaled=F, # influence could be expressed in other units
                      scale.by=median # e.g. percent of median influence etc
  ){
    diff = apply(M_i,margin,function(x){sum(abs(x-tariff))})
    if(scaled){
      diff= diff /scale.by(diff)
    }
    diff
  }
  
  