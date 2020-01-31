
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
  
  # retrieved diserpsion measures from 2 social value sets
  get_dispersion = function(CAM,MUAP){
    df = data.frame(mean = c(mean(CAM),mean(MUAP)),
                    sd = c(sd(CAM),sd(MUAP)),
                    min = c(min(CAM),min(MUAP)),
                    max = c(max(CAM)[1],max(MUAP)[1]),
                    range = c(max(CAM)[1]-min(CAM)[1],max(MUAP)[1]-min(MUAP)[1]),
                    max_min_ratio = c(max(CAM)[1]/min(CAM)[1],max(MUAP)[1]/min(MUAP)[1]),
                    q25 = c(quantile(CAM,0.25),quantile(MUAP,0.25)),
                    median = c(median(CAM),median(MUAP)),
                    q75 = c(quantile(CAM,0.75),quantile(MUAP,0.75)),
                    iqr = c(IQR(CAM),IQR(MUAP))
                    
    )
    df = data.frame(t(df))
    names(df) = c("CAM","MUAP")
    df = cbind(df,ratio=df[,"MUAP"]/df[,"CAM"])
    df$ratio = ifelse(c(F,T,F,F,T,T,F,F,F,T),df$ratio,NA)
    df = data.frame(apply(df,2,function(x){formatC(x,digits = 3,format = "f")}),stringsAsFactors = F)
    df$ratio[grepl("NA",df$ratio)]  = ""
    df
  }
  
  # extracts data from lm models
  make_tbl = function(model,get.rownames=F,give.colname=NULL){
    temp = formatC(cbind(
      model$coefficients,confint(model)[,1],confint(model)[,2]),digits=2,format="f")  
    temp = paste(temp[,1]," (",temp[,2],"; ",temp[,3],")",sep="")
    r2 = formatC(summary(model)$r.squared,digits=2,format="f")
    n = summary(model)$df[1] +summary(model)$df[2]
    temp = c(temp,r2,n)
    if(get.rownames){
      temp = cbind(c(names(model$coefficients),"R2","Observations"),
                   temp)
      temp[,1] = gsub("TRUE","",temp[,1])
    } else{
      temp = cbind(temp)
      colnames(temp) = give.colname
    }
    temp
  }
  
  # extracts data from lm models
  make_mtbl = function(m){ 
    temp = data.frame(m$coefficients,confint(m))
    temp = apply(temp,2,function(x){formatC(x,digits=2,format="f")})
    temp = apply(temp,1,function(x){paste(x[1]," (",x[2],"; ",x[3],")",sep="")})
    temp = data.frame(cbind(temp),stringsAsFactors = F)
    names(temp) = ""
    if(length(temp[,1])==2){
      temp = rbind(temp,c(""),c(""))
    }
    if(length(temp[,1])==3){
      temp = rbind(temp,c(""))
    }
    
    temp = rbind(temp,"R$^2$" = formatC(summary(m)$r.squared,digits=2,format="f"),"N"=summary(m)$df[1]+summary(m)$df[2])
    
    return(temp)
  }  

  # computes individuals' influences
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
  