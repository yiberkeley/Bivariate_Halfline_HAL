
#Poisson cv fit helper function
poisson_cv_fit_func_cvxr<-function(poisson_data_N3,basis_list_N3,nfolds,penalty.factor_N3,undersmooth=F,l1_norm_start=F,uniform_ratio=F,
                                   offset_initial=NULL,length_clever=NULL,coef_initial=NULL,offset_relax=NULL){
  #Implemented the relaxed targeting
  if(is.null(offset_relax) | sum(offset_relax==F)>=1){
      #CV part
      if(l1_norm_start==F & is.logical(l1_norm_start)){
        low<-1
        up<-10
        l1_norm_seq<-seq(low,up,length.out=5)
        print(l1_norm_start)
      }else{
        low<-max(l1_norm_start-1,0)
        up<-l1_norm_start+1
        l1_norm_seq<-seq(low,up,length.out=5)
        print(l1_norm_start)
      }
      stop_cv<-FALSE
      zoom_in_times<-3
      zoom_in_middle<-0
      while(stop_cv==FALSE & zoom_in_times <6){
        print("check point: prepare to get cv poisson fit")

        print("check point test 1" %+% length_clever)
        fitness<-cv_fit_cvxr(poisson_data_N3,l1_norm_seq,nfolds,penalty.factor_N3,uniform_ratio,
                             position_basis=NULL,additional_offset=NULL,offset_initial,length_clever,coef_initial)

        print("check point: cv poisson fit finished")
        print(fitness$l1_norm)
        l1_norm_cv<-as.numeric(fitness$l1_norm[which(fitness$cv_loss==min(fitness$cv_loss))])

        if ( sum((unique(l1_norm_cv)==sort(fitness$l1_norm)[1]) |
                 (unique(l1_norm_cv)==sort(fitness$l1_norm)[1]))>=1 ){
          low_temp<-min(l1_norm_cv)
          up<-low_temp+1/zoom_in_times
          low<-max(low_temp-1/zoom_in_times,0)
          l1_norm_seq<-seq(up,low,length.out=5)
          zoom_in_times<-zoom_in_times+1
          print("cv zoom in low")
        }else if (sum((unique(l1_norm_cv)==sort(fitness$l1_norm,decreasing = T)[1]) |
                      (unique(l1_norm_cv)==sort(fitness$l1_norm,decreasing = T)[1]))>=1 ){

          up_temp<-max(l1_norm_cv)
          up<-up_temp+1/zoom_in_times
          low<-max(up_temp-1/zoom_in_times,0)
          l1_norm_seq<-seq(low,up,length.out=5)
          zoom_in_times<-zoom_in_times+1
          print("cv zoom in up")
        }else if (min(l1_norm_cv)==0){
          stop_cv<-TRUE
          print("cv finished with intercept selected only")
        }else if (zoom_in_middle<=1 & l1_norm_start==F & is.logical(l1_norm_start)){ #zoom_in_middle<=2
          middle_temp<-min(l1_norm_cv)
          up<-middle_temp+1/zoom_in_times
          low<-max(middle_temp-1/zoom_in_times,0)
          l1_norm_seq<-seq(up,low,length.out=5)
          zoom_in_times<-zoom_in_times+1
          zoom_in_middle<-zoom_in_middle+1
          print("cv zoom in middle")
        }else if (zoom_in_middle<=0 & !(l1_norm_start==F & is.logical(l1_norm_start)) ){ #zoom_in_middle<=1
          middle_temp<-min(l1_norm_cv)
          up<-middle_temp+1/zoom_in_times
          low<-max(middle_temp-1/zoom_in_times,0)
          l1_norm_seq<-seq(up,low,length.out=5)
          zoom_in_times<-zoom_in_times+1
          zoom_in_middle<-zoom_in_middle+1
          print("cv zoom in middle")
        }else{
          stop_cv<-TRUE
          print("cv finished")
        }
      }

      l1_norm_chosen_seq=min((fitness$l1_norm[which(fitness$cv_loss==min(fitness$cv_loss))])) #only cv
      cv_selected_cv_risk_N3<-min(fitness$cv_loss)
      cv_selected_cv_number_coef_N3<-unique(fitness$cv_num_coef[which(fitness$cv_loss==min(fitness$cv_loss))])


      offset<-log(poisson_data_N3$timeIntervalLength)
      weights<- poisson_data_N3$weights
      y<-poisson_data_N3$jump
      x_basis_poisson<-as.matrix(poisson_data_N3[,-c(1:9)])
      time_vec<-as.vector(poisson_data_N3$t)

      if(is.null(offset_initial)){
        fit_N3<-cvrx_poisson_l1(l1_norm_chosen_seq,x_basis_poisson,offset,weights,y,penalty.factor_N3,time_vec,uniform_ratio)
      }else{
        fit_N3<-cvrx_poisson_l1_offset(l1_norm_chosen_seq,x_basis_poisson,offset,weights,y,penalty.factor_N3,time_vec,uniform_ratio,length_clever,coef_initial)
      }
  }else{

    print("#####enter the right offset region with offset relax as TRUE")
    l1_norm_chosen_seq=100
    cv_selected_cv_risk_N3<-(-1)
    cv_selected_cv_number_coef_N3<-(-1)

    offset<-log(poisson_data_N3$timeIntervalLength)
    weights<- poisson_data_N3$weights
    y<-poisson_data_N3$jump
    x_basis_poisson<-as.matrix(poisson_data_N3[,-c(1:9)])
    time_vec<-as.vector(poisson_data_N3$t)

    fit_N3<-cvrx_poisson_l1_offset(100,x_basis_poisson,offset,weights,y,penalty.factor_N3,time_vec,uniform_ratio,length_clever,coef_initial)
  }



  beta<-fit_N3[[1]][[1]]
  beta_0<-fit_N3[[1]][[2]]
  cv_selected_fit_risk_N3<-fit_N3[[1]][[3]]

  coef_N3_initial_list<-list()
  basis_list_N3_select_list<-list()
  penalty_N3_select_list<-list()
  for(i in 1:length(l1_norm_chosen_seq)){
    coef_N3_initial_list[[i]]<-c(beta_0,beta[abs(beta)>10^-4])
    basis_list_N3_select_list[[i]]<-basis_list_N3[abs(beta)>10^-4]
    penalty_N3_select_list[[i]]<-penalty.factor_N3[abs(beta)>10^-4]
  }


  #Undersmoothing part
  if(undersmooth==T){

      low<-l1_norm_chosen_seq
      up<-l1_norm_chosen_seq+2
      l1_norm_seq<-seq(low,up,length.out=5)
      print(l1_norm_start)

      stop_undersmooth<-FALSE
      zoom_in_times<-1

      while(stop_undersmooth==FALSE & zoom_in_times <6){ #10
        print("check point: prepare to get undersmooth cvxr fit")
        fitness<-cv_fit_cvxr(poisson_data_N3,l1_norm_seq,nfolds,penalty.factor_N3,uniform_ratio,
                             position_basis=NULL,additional_offset=NULL,offset_initial,length_clever,coef_initial)
        print("check point: undersmooth poisson fit finished")
        if ( sum(fitness$cv_num_coef>cv_selected_cv_number_coef_N3*3)==0){ #1.5 #| sum(fitness$cv_loss>2*cv_selected_cv_risk_N3)==0
          low<-up
          up<-up+2
          l1_norm_seq<-seq(low,up,length.out=5)
          zoom_in_times<-zoom_in_times+1
          print("undersmooth zoom in up")
        }else{
          stop_undersmooth<-TRUE
          print("undersmooth finished")
        }
      }

      lambda_chosen_undersmooth=max(exp(fitness$l1_norm[which( fitness$cv_num_coef>cv_selected_cv_number_coef_N3*1.5)]))

    l1_norm_chosen_undersmooth=fitness$l1_norm[min(which(fitness$cv_num_coef>cv_selected_cv_number_coef_N3*1.5),length(fitness))]



    offset<-log(poisson_data_N3$timeIntervalLength)
    weights<- poisson_data_N3$weights
    y<-poisson_data_N3$jump
    x_basis_poisson<-as.matrix(poisson_data_N3[,-c(1:9)])
    time_vec<-as.vector(poisson_data_N3$t)


    if(is.null(offset_initial)){
      fit_N3_undersmooth<-cvrx_poisson_l1(l1_norm_chosen_undersmooth,x_basis_poisson,offset,weights,y,penalty.factor_N3,time_vec,uniform_ratio)
    }else{
      fit_N3_undersmooth<-cvrx_poisson_l1_offset(l1_norm_chosen_undersmooth,x_basis_poisson,offset,weights,y,penalty.factor_N3,time_vec,uniform_ratio,length_clever,coef_initial)
    }


    beta_undersmooth<-fit_N3_undersmooth[[1]][[1]]
    beta_0_undersmooth<-fit_N3_undersmooth[[1]][[2]]

    coef_N3_initial_list_undersmooth<-list()
    basis_list_N3_select_list_undersmooth<-list()
    penalty_N3_select_list_undersmooth<-list()
    for(i in 1:length(l1_norm_chosen_seq)){
      coef_N3_initial_list_undersmooth[[i]]<-c(beta_0_undersmooth,beta_undersmooth[abs(beta_undersmooth)>10^-4])
      basis_list_N3_select_list_undersmooth[[i]]<-basis_list_N3[abs(beta_undersmooth)>10^-4]
      penalty_N3_select_list_undersmooth[[i]]<-penalty.factor_N3[abs(beta_undersmooth)>10^-4]
    }

  }


  result_list<-list()
  result_list[[1]]<-cv_selected_cv_risk_N3
  result_list[[2]]<-cv_selected_cv_number_coef_N3
  result_list[[3]]<-cv_selected_fit_risk_N3
  result_list[[4]]<-coef_N3_initial_list
  result_list[[5]]<-basis_list_N3_select_list
  result_list[[6]]<-penalty_N3_select_list
  if(undersmooth==T){
    result_list[[7]]<-coef_N3_initial_list_undersmooth
    result_list[[8]]<-basis_list_N3_select_list_undersmooth
    result_list[[9]]<-penalty_N3_select_list_undersmooth
  }
  return(result_list)
}
