#in build, install
library(devtools)
library(hal9001)
load_all()
library(glmnet)
library(dplyr)
library(CVXR)

`%+%` <- function(a, b) paste0(a, b)

library(glmnet)
library(Rlab)

library(foreach)
library(doParallel)

library(future.apply)
library(fMultivar)

library(doFuture)

run_sim<-function(n,
                 seed,
                 df,
                 nfolds,
                 undersmooth_type="cv_all",
                 basis_option="0_order",
                 masking_option="no_masking",
                 weight_option="no_tail_shrinkage",
                 penalty_type="usual",
                 undersmooth=F,
                 censoring=F,
                 perc_sample=0.05,
                 perc_sample_basis=0.05,
                 cvxr=T,
                 uniform_ratio=F){


  data_full<-df
  data_obs<-data_full%>%select(tt1,delt1,tt2,delt2)



  ##Set up hyper parameters
  #Create a fold id for future cross evaluation
  nfolds<-nfolds
  foldid = sample(rep(seq(nfolds), length = n))
  #Create the subsampling percentage for monitoring time
  perc_sample=perc_sample
  #Create the subsampling percentage for the basis
  perc_sample_basis=perc_sample_basis



  ##Set up the monitoring time for repeated data
  #Monitoring time, starting from the smallest observed event time - sd(observed event times)
  monitoring_time<-sort(unique(c(data_obs$tt1,data_obs$tt2)))
  #monitoring_time<-c(min(monitoring_time)-sd(monitoring_time),monitoring_time)
  monitoring_time<-unique(c(0,monitoring_time))




  ##Created repeated data
  repeated_data_collection<-repeated_data_func(data_obs,foldid,perc_sample,monitoring_time,censor_data_include=T,masking_option,weight_option)
  repeated_data_N1_2<-repeated_data_collection[[1]]
  repeated_data_N2_2<-repeated_data_collection[[2]]
  if(censoring){
    repeated_data_A1_2<-repeated_data_collection[[3]]
    repeated_data_A2_2<-repeated_data_collection[[4]]
  }



  ##Prepare the poisson data and do the poisson fit
  if(basis_option=="0_order"){
    #Augment the Repeated data with basis
    #Basis list set up for N1,A1
    basis_list_N1 <- enumerate_basis(repeated_data_N1_2%>%select(t,tt2_less_t,tt2_less_t_delt2,tt2_less_t_tt2),max_degree = 2,smoothness_orders = 0)
    if(censoring){
      basis_list_A1 <- enumerate_basis(repeated_data_A1_2%>%select(t,tt2_less_t,tt2_less_t_delt2,tt2_less_t_tt2),max_degree = 2,smoothness_orders = 0)
    }
    #Subsample basis list
    basis_list_N1<-basis_list_N1[sample(1:length(basis_list_N1),floor(perc_sample_basis*length(basis_list_N1)),replace=FALSE)]
    if(censoring){
      basis_list_A1<-basis_list_A1[sample(1:length(basis_list_A1),floor(perc_sample_basis*length(basis_list_A1)),replace=FALSE)]
    }

    #Basis list set up for N2
    basis_list_N2 <- enumerate_basis(repeated_data_N2_2%>%select(t,tt1_less_t,tt1_less_t_delt1,tt1_less_t_tt1),max_degree = 2,smoothness_orders = 0)
    if(censoring){
      basis_list_A2 <- enumerate_basis(repeated_data_A2_2%>%select(t,tt1_less_t,tt1_less_t_delt1,tt1_less_t_tt1),max_degree = 2,smoothness_orders = 0)
    }
    #Subsample basis list
    basis_list_N2<-basis_list_N2[sample(1:length(basis_list_N2),floor(perc_sample_basis*length(basis_list_N2)),replace=FALSE)]
    if(censoring){
      basis_list_A2<-basis_list_A2[sample(1:length(basis_list_A2),floor(perc_sample_basis*length(basis_list_A2)),replace=FALSE)]
    }
  }else if(basis_option=="1_order"){
    #Augment the Repeated data with basis
    #Basis list set up for N1,A1
    basis_list_N1 <- enumerate_basis(repeated_data_N1_2%>%select(t,tt2_less_t,tt2_less_t_delt2,tt2_less_t_tt2),max_degree = 2,smoothness_orders = 1)
    if(censoring){
      basis_list_A1 <- enumerate_basis(repeated_data_A1_2%>%select(t,tt2_less_t,tt2_less_t_delt2,tt2_less_t_tt2),max_degree = 2,smoothness_orders = 1)
    }
    #Subsample basis list
    basis_list_N1<-basis_list_N1[sample(1:length(basis_list_N1),floor(perc_sample_basis*length(basis_list_N1)),replace=FALSE)]
    if(censoring){
      basis_list_A1<-basis_list_A1[sample(1:length(basis_list_A1),floor(perc_sample_basis*length(basis_list_A1)),replace=FALSE)]
    }

    #Basis list set up for N2
    basis_list_N2 <- enumerate_basis(repeated_data_N2_2%>%select(t,tt1_less_t,tt1_less_t_delt1,tt1_less_t_tt1),max_degree = 2,smoothness_orders = 1)
    if(censoring){
      basis_list_A2 <- enumerate_basis(repeated_data_A2_2%>%select(t,tt1_less_t,tt1_less_t_delt1,tt1_less_t_tt1),max_degree = 2,smoothness_orders = 1)
    }
    #Subsample basis list
    basis_list_N2<-basis_list_N2[sample(1:length(basis_list_N2),floor(perc_sample_basis*length(basis_list_N2)),replace=FALSE)]
    if(censoring){
     basis_list_A2<-basis_list_A2[sample(1:length(basis_list_A2),floor(perc_sample_basis*length(basis_list_A2)),replace=FALSE)]
    }
  }


  penalty_rep_func<-function(basis_list_N3){
    position_basis_N3_1<-lapply(basis_list_N3, `[[`, 1)
    position_basis_N3_1<-lapply(position_basis_N3_1,function(x){sum(unlist(x)==4)>=1 | sum(unlist(x)==1)>=1})
    position_basis_N3_1<-which(unlist(position_basis_N3_1)!=TRUE)

    penalty_temp<-lapply(basis_list_N3, `[[`, 2)
    penalty_temp<-lapply(penalty_temp,function(x){
      if(length(x)==2 & sum(x!=0 &x!=1)!=0){
        1/(max(x[which(x!=0 &x!=1)])+0.000001)
      }else{
        1/(x[1]+0.000001)
      }
    })

    penalty_temp<-unlist(penalty_temp)
    penalty_temp[position_basis_N3_1]<-2
    penalty_temp<-(penalty_temp/sum(penalty_temp))*length(penalty_temp)
    return(penalty_temp)
  }


  ##Penalty.factor set up for poisson fit
  #Usual scenario
  if(penalty_type=="usual"){
    penalty.factor_N1<-rep(1,length(basis_list_N1))
    if(censoring){
      penalty.factor_A1<-rep(1,length(basis_list_A1))
    }
    penalty.factor_N2<-rep(1,length(basis_list_N2))
    if(censoring){
      penalty.factor_A2<-rep(1,length(basis_list_A2))
    }
  }

  if(penalty_type=="rep"){
    penalty.factor_N1<-penalty_rep_func(basis_list_N1)
    if(censoring){
      penalty.factor_A1<-penalty_rep_func(basis_list_A1)
    }
    penalty.factor_N2<-penalty_rep_func(basis_list_N2)
    if(censoring){
      penalty.factor_A2<-penalty_rep_func(basis_list_A2)
    }
  }
  #Flex t scenario
  if(penalty_type=="flex_t"){
    #N1 jumping process
    reult_temp<-flex_t_func(data_obs,basis_list_N1,list_type="N1",basis_option)
    basis_list_N1<-reult_temp[[1]]
    penalty.factor_N1<-reult_temp[[2]]

    if(censoring){
      #A1 jumping process
      reult_temp<-flex_t_func(data_obs,basis_list_A1,list_type="A1",basis_option)
      basis_list_A1<-reult_temp[[1]]
      penalty.factor_A1<-reult_temp[[2]]
    }

    #N2 jumping process
    reult_temp<-flex_t_func(data_obs,basis_list_N2,list_type="N2",basis_option)
    basis_list_N2<-reult_temp[[1]]
    penalty.factor_N2<-reult_temp[[2]]

    if(censoring){
      #A2 jumping process
      reult_temp<-flex_t_func(data_obs,basis_list_A2,list_type="A2",basis_option)
      basis_list_A2<-reult_temp[[1]]
      penalty.factor_A2<-reult_temp[[2]]
    }
  }

  #Flex tail region scenerio
  if(penalty_type=="flex_tail"){
    #N1 jumping process
    reult_temp<-flex_tail_func(data_obs,basis_list_N1,list_type="N1",basis_option)
    basis_list_N1<-reult_temp[[1]]
    penalty.factor_N1<-reult_temp[[2]]

    if(censoring){
      #A1 jumping process
      reult_temp<-flex_tail_func(data_obs,basis_list_A1,list_type="A1",basis_option)
      basis_list_A1<-reult_temp[[1]]
      penalty.factor_A1<-reult_temp[[2]]
    }

    #N2 jumping process
    reult_temp<-flex_tail_func(data_obs,basis_list_N2,list_type="N2",basis_option)
    basis_list_N2<-reult_temp[[1]]
    penalty.factor_N2<-reult_temp[[2]]

    if(censoring){
      #A2 jumping process
      reult_temp<-flex_tail_func(data_obs,basis_list_A2,list_type="A2",basis_option)
      basis_list_A2<-reult_temp[[1]]
      penalty.factor_A2<-reult_temp[[2]]
    }
  }




  ##Now,naugment the Repeated data with basis
  #Poisson data N1,A1
  basis_expansion<-make_design_matrix(as.matrix(repeated_data_N1_2%>%select(t,tt2_less_t,tt2_less_t_delt2,tt2_less_t_tt2)),basis_list_N1)
  basic_info<-repeated_data_N1_2
  poisson_data_N1<-cbind(basic_info,as.matrix(basis_expansion))

  if(censoring){
    basis_expansion<-make_design_matrix(as.matrix(repeated_data_A1_2%>%select(t,tt2_less_t,tt2_less_t_delt2,tt2_less_t_tt2)),basis_list_A1)
    basic_info<-repeated_data_A1_2
    poisson_data_A1<-cbind(basic_info,as.matrix(basis_expansion))
  }
  #Poisson data N2,A2
  basis_expansion<-make_design_matrix(as.matrix(repeated_data_N2_2%>%select(t,tt1_less_t,tt1_less_t_delt1,tt1_less_t_tt1)),basis_list_N2)
  basic_info<-repeated_data_N2_2
  poisson_data_N2<-cbind(basic_info,as.matrix(basis_expansion))

  if(censoring){
    basis_expansion<-make_design_matrix(as.matrix(repeated_data_A2_2%>%select(t,tt1_less_t,tt1_less_t_delt1,tt1_less_t_tt1)),basis_list_A2)
    basic_info<-repeated_data_A2_2
    poisson_data_A2<-cbind(basic_info,as.matrix(basis_expansion))
  }


  poisson_cv_fit_func_all<-function(cvxr,poisson_data_N3,basis_list_N3,nfolds,penalty.factor_N3,undersmooth,l1_norm_start=F,uniform_ratio=F){
    if(!cvxr){
      result<-poisson_cv_fit_func(poisson_data_N3,basis_list_N3,nfolds,penalty.factor_N3,undersmooth,l1_norm_start,uniform_ratio)
    }else{
      result<-poisson_cv_fit_func_cvxr(poisson_data_N3,basis_list_N3,nfolds,penalty.factor_N3,undersmooth,l1_norm_start,uniform_ratio)
    }
    return(result)
  }


  ##Continue poisson fit

  if(undersmooth_type=="cv_all"){
      #N1 jumping process
      result_temp<-poisson_cv_fit_func_all(cvxr,poisson_data_N1,basis_list_N1,nfolds,penalty.factor_N1,undersmooth,l1_norm_start=F,uniform_ratio)
      cv_selected_cv_risk_N1 <- result_temp[[1]]
      cv_selected_cv_number_coef_N1 <- result_temp[[2]]
      cv_selected_fit_risk_N1 <- result_temp[[3]]
      coef_N1_initial_list <- result_temp[[4]]
      basis_list_N1_select_list <- result_temp[[5]]
      penalty_N1_select_list <- result_temp[[6]]
      print("N1 fit finished")
      if(undersmooth==T){
        coef_N1_initial_list_undersmooth <- result_temp[[7]]
        basis_list_N1_select_list_undersmooth <- result_temp[[8]]
        penalty_N1_select_list_undersmooth <- result_temp[[9]]
        print("N1 undersmooth fit finished")
      }

      if(censoring){
        #A1 jumping process
        result_temp<-poisson_cv_fit_func_all(cvxr,poisson_data_A1,basis_list_A1,nfolds,penalty.factor_A1,undersmooth,l1_norm_start=F,uniform_ratio)
        cv_selected_cv_risk_A1 <- result_temp[[1]]
        cv_selected_cv_number_coef_A1 <- result_temp[[2]]
        cv_selected_fit_risk_A1 <- result_temp[[3]]
        coef_A1_initial_list <- result_temp[[4]]
        basis_list_A1_select_list <- result_temp[[5]]
        penalty_A1_select_list <- result_temp[[6]]
        print("A1 fit finished")
        if(undersmooth==T){
          coef_A1_initial_list_undersmooth <- result_temp[[7]]
          basis_list_A1_select_list_undersmooth <- result_temp[[8]]
          penalty_A1_select_list_undersmooth <- result_temp[[9]]
          print("A1 undersmooth fit finished")
        }
      }

      #N2 jumping process
      result_temp<-poisson_cv_fit_func_all(cvxr,poisson_data_N2,basis_list_N2,nfolds,penalty.factor_N2,undersmooth,l1_norm_start=F,uniform_ratio)
      cv_selected_cv_risk_N2 <- result_temp[[1]]
      cv_selected_cv_number_coef_N2 <- result_temp[[2]]
      cv_selected_fit_risk_N2 <- result_temp[[3]]
      coef_N2_initial_list <- result_temp[[4]]
      basis_list_N2_select_list <- result_temp[[5]]
      penalty_N2_select_list <- result_temp[[6]]
      print("N2 fit finished")
      if(undersmooth==T){
        coef_N2_initial_list_undersmooth <- result_temp[[7]]
        basis_list_N2_select_list_undersmooth <- result_temp[[8]]
        penalty_N2_select_list_undersmooth <- result_temp[[9]]
        print("N2 undersmooth fit finished")
      }


      if(censoring){
        #A2 jumping process
        result_temp<-poisson_cv_fit_func_all(cvxr,poisson_data_A2,basis_list_A2,nfolds,penalty.factor_A2,undersmooth,l1_norm_start=F,uniform_ratio)
        cv_selected_cv_risk_A2 <- result_temp[[1]]
        cv_selected_cv_number_coef_A2 <- result_temp[[2]]
        cv_selected_fit_risk_A2 <- result_temp[[3]]
        coef_A2_initial_list <- result_temp[[4]]
        basis_list_A2_select_list <- result_temp[[5]]
        penalty_A2_select_list <- result_temp[[6]]
        print("A2 fit finished")
        if(undersmooth==T){
          coef_A2_initial_list_undersmooth <- result_temp[[7]]
          basis_list_A2_select_list_undersmooth <- result_temp[[8]]
          penalty_A2_select_list_undersmooth <- result_temp[[9]]
          print("A2 undersmooth fit finished")
        }
      }
  }


      coef_N1_initial = coef_N1_initial_list[[1]]
      basis_list_N1_select = basis_list_N1_select_list[[1]]
      coef_N2_initial = coef_N2_initial_list[[1]]
      basis_list_N2_select = basis_list_N2_select_list[[1]]
      coef_A1_initial = ifelse(censoring,coef_A1_initial_list[[1]],-1)
      basis_list_A1_select = ifelse(censoring,basis_list_A1_select_list[[1]],-1)
      coef_A2_initial =ifelse(censoring, coef_A2_initial_list[[1]],-1)
      basis_list_A2_select = ifelse(censoring,basis_list_A2_select_list[[1]],-1)

      result<-list()
      result[[1]]<-coef_N1_initial
      result[[2]]<-basis_list_N1_select
      result[[3]]<-coef_N2_initial
      result[[4]]<-basis_list_N2_select
      result[[5]]<-coef_A1_initial
      result[[6]]<-basis_list_A1_select
      result[[7]]<-coef_A2_initial
      result[[8]]<-basis_list_A2_select
    return(result)
}



