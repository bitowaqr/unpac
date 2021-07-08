### Simple example to motivate the UNPAc method
  # utilities for health states full health, A, B, AB
  Alice = c(1,0.966666666,0.9333333,0.9)
  Bob = c(1,-0.333333333,0.333333333,-1)
  ase_df = (t(cbind(Alice,Bob)))
  colnames(ase_df) = c("h_full","h_a","h_b","h_ab")
  # standard social value set
  ase_standard = colMeans(ase_df)
  # normalised utilities
  ase_norm_utils = data.frame(t(apply(ase_df,1,function(x){(x-min(x))/(max(x)-min(x))})))
  # relative social value set
  ase_rel_unpac = colMeans(ase_norm_utils)
  # utility ranges
  ase_util_ranges = unlist(apply(ase_df,1,function(x){max(x)-min(x)}))
  # scaling factor
  ase_scaling_f = mean(ase_util_ranges)
  # unpac social value set
  ase_unpac = ase_rel_unpac * ase_scaling_f + (1-ase_scaling_f)
  
  p_ase = 
    ggplot() +
    geom_hline(yintercept = 0,size=0.3) +
    geom_point(aes(x=1:4,y=ase_df[1,],col="Alice"),size=2) +
    geom_path(aes(x=1:4,y=ase_df[1,],col="Alice",linetype="Alice"),size=1) +
    geom_point(aes(x=1:4,y=ase_df[2,],col="Bob"),size=2) +
    geom_path(aes(x=1:4,y=ase_df[2,],col="Bob",linetype="Bob"),size=1) +
    geom_point(aes(x=1:4,y=ase_standard,col="standard"),size=2) +
    geom_path(aes(x=1:4,y=ase_standard,linetype="standard",col="standard"),size=1.1) +
    geom_point(aes(x=1:4,y=ase_unpac,col="unpac"),size=2) +
    geom_path(aes(x=1:4,y=ase_unpac,linetype="unpac",col="unpac"),size=1.1) +
    scale_x_continuous(labels = colnames(ase_df), name="Health state") +
    ylab("Utility value")  +
    theme_minimal() + 
    scale_linetype_discrete(guide = FALSE) +
    scale_color_brewer(palette = "Dark2",name = "") 
  
  ggsave(filename = "./output/p_ase.png",plot=p_ase,height = 5,width=7)
