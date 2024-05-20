#Once get the conditional intensity estimated, we want to get the canonical gradient and the estimated survival curve.

Rejection_Sample_CdfBased<-function(coef_N1_initial,basis_list_N1_select,
                                        coef_N2_initial,basis_list_N2_select,
                                        coef_A1_initial,basis_list_A1_select,
                                        coef_A2_initial,basis_list_A2_select,
                                        censoring,
                                        half_line_position,
                                        half_line_type,
                                        starting_of_half_line,
                                        sample_size
){

  t_cut=half_line_position
  t_start=starting_of_half_line
  coef_N1<-coef_N1_initial
  coef_N2<-coef_N2_initial

  #Extact the time points where the fitted conditional intensity is possible to change
  #N1 jumping process
  col_index_temp<-unlist(lapply(basis_list_N1_select, `[[`, 1))
  val_temp<-unlist(lapply(basis_list_N1_select, `[[`, 2))
  N1_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
  #N2 jumping process
  col_index_temp<-unlist(lapply(basis_list_N2_select, `[[`, 1))
  val_temp<-unlist(lapply(basis_list_N2_select, `[[`, 2))
  N2_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
  if(censoring){
    #A1 jumping process
    col_index_temp<-unlist(lapply(basis_list_A1_select, `[[`, 1))
    val_temp<-unlist(lapply(basis_list_A1_select, `[[`, 2))
    A1_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
    #A2 jumping process
    col_index_temp<-unlist(lapply(basis_list_A2_select, `[[`, 1))
    val_temp<-unlist(lapply(basis_list_A2_select, `[[`, 2))
    A2_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
  }else{
    A1_time_grid<-NULL
    A2_time_grid<-NULL
  }

  #Get the cdf

  single_censored_prob<-function(x){
    prob_overall<-1
    if(half_line_type=="N1"){
      if(x<=t_cut){
        prob_overall<-prob_overall*product_integral_func(start=0,end=x,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
        prob_overall<-prob_overall*product_integral_func(start=0,end=x,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
        prob_overall<-prob_overall*product_integral_func(start=x,end=t_cut,1,0,x,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
        prob_overall<-prob_overall*hazard_func(t_cut,1,0,x,basis_list_N2_select,coef_N2,type="N2")

      }else{
        prob_overall<-prob_overall*product_integral_func(start=0,end=x,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
        prob_overall<-prob_overall*product_integral_func(start=0,end=x,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
        prob_overall<-prob_overall*hazard_func(t_cut,0,0,0,basis_list_N2_select,coef_N2,type="N2")
        prob_overall<-prob_overall*product_integral_func(start=t_cut,end=x,1,1,t_cut,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
      }
    }else{
      if(x<=t_cut){
        prob_overall<-prob_overall*product_integral_func(start=0,end=x,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
        prob_overall<-prob_overall*product_integral_func(start=0,end=x,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
        prob_overall<-prob_overall*product_integral_func(start=x,end=t_cut,1,0,x,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
        prob_overall<-prob_overall*hazard_func(t_cut,1,0,x,basis_list_N1_select,coef_N1,type="N1")

      }else{
        prob_overall<-prob_overall*product_integral_func(start=0,end=x,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
        prob_overall<-prob_overall*product_integral_func(start=0,end=x,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
        prob_overall<-prob_overall*hazard_func(t_cut,0,0,0,basis_list_N1_select,coef_N1,type="N1")
        prob_overall<-prob_overall*product_integral_func(start=t_cut,end=x,1,1,t_cut,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
      }
    }
  }

  #Use rejection sampling technique
  grid_scale_density<-0.01
  grid_temp<- seq(t_start,1,grid_scale_density)

  single_censored_prob_values<-unlist(lapply(grid_temp,function(x)single_censored_prob(x)))

  cdf_values<-(1-single_censored_prob_values/single_censored_prob_values[1])
  cdf_values[length(cdf_values)]<-1

  cdf_sampled<-  runif(sample_size,0,1)

  inverse_cdf<-function(x){
    y_min<-max(cdf_values[cdf_values<=x])
    x_min<-max(grid_temp[cdf_values<=x])

    y_max<-min(cdf_values[cdf_values>x])
    x_max<-min(grid_temp[cdf_values>x])

    ratio<-(x-y_min)/(y_max-y_min)
    draw<-x_min+ratio*(x_max-x_min)
    weight<-(y_max-y_min)/(x_max-x_min)
    return(c(draw,weight))
  }

  final_sample<-NULL
  final_sample_weight<-NULL
  for(i in 1:sample_size){
    result_temp<-inverse_cdf(cdf_sampled[i])
    final_sample<-c(final_sample,result_temp[1])
    final_sample_weight<-c(final_sample_weight,result_temp[2])
  }

  result_inner<-list()
  result_inner[[1]]<-final_sample
  result_inner[[2]]<-final_sample_weight
  return(result_inner)
}


