source("utils.R")


#n=500, perc_sample=0.1,perc_sample_basis=0.1
#n=1000, monitoring times: perc_sample=0.02 & full basis subsampling: perc_sample_basis=0.02
#n=2000, perc_sample=0.01,perc_sample_basis=0.01
#seed_vals<-c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)
seed_vals<-c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)
seed<-seed_vals[1]
  n=100
  res <- run_sim(n = n,
                 seed = seed,
                 df=df,
                 nfolds=3,
                 undersmooth_type="cv_all",
                 basis_option="0_order",
                 masking_option="no_masking", #0.7
                 weight_option="no_tail_shrinkage",
                 penalty_type="usual",
                 undersmooth=F,
                 censoring=F,
                 perc_sample=0.05,
                 perc_sample_basis=0.05,
                 cvxr=T,
                 uniform_ratio=F)

  Sample<-Rejection_Sample_DensityBased(coef_N1_initial = res[[1]],
                           basis_list_N1_select = res[[2]],
                           coef_N2_initial = res[[3]],
                           basis_list_N2_select = res[[4]],
                           coef_A1_initial = res[[5]],
                           basis_list_A1_select = res[[6]],
                           coef_A2_initial =res[[7]],
                           basis_list_A2_select = res[[8]],
                           censoring=F,
                           half_line_position=0.7, #where the observed time is
                           half_line_type="N2", #which the censored event is
                           starting_of_half_line=0.2, # where the censored event happen
                           sample_size=5
  )

  Sample2<-Rejection_Sample_CdfBased(coef_N1_initial = res[[1]],
                                        basis_list_N1_select = res[[2]],
                                        coef_N2_initial = res[[3]],
                                        basis_list_N2_select = res[[4]],
                                        coef_A1_initial = res[[5]],
                                        basis_list_A1_select = res[[6]],
                                        coef_A2_initial =res[[7]],
                                        basis_list_A2_select = res[[8]],
                                        censoring=F,
                                        half_line_position=0.7, #where the observed time is
                                        half_line_type="N2", #which the censored event is
                                        starting_of_half_line=0.2, # where the censored event happen
                                        sample_size=5
  )

