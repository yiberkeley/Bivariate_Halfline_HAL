##Some helpful functions
product_integral_func<-function(start,end,first,second,third,N3_time_grid,basis_list_N3_select,coef_N3,type){
  prob<-1
  if(start<end){
    if(type=="A1" |type=="N1"){
      grid_1_N3<-sort(c(N3_time_grid[N3_time_grid<end & N3_time_grid>start],start,end))
      for(i in 1:(length(grid_1_N3)-1)){
        temp<-data.frame(t=(grid_1_N3[i]+grid_1_N3[i+1])/2,tt2_less_t=first,tt2_less_t_delt2=second,tt2_less_t_tt2=third)
        temp2<-hal9001::make_design_matrix(as.matrix(temp),basis_list_N3_select)
        temp2<-c(1,as.vector(temp2))
        log_hazard_midpoint_N3<-sum(coef_N3*temp2)
        prob_N3_interval<-exp(-(grid_1_N3[i+1]-grid_1_N3[i])*exp(log_hazard_midpoint_N3))
        prob<-prob*prob_N3_interval
      }
    }else{
      grid_1_N3<-sort(c(N3_time_grid[N3_time_grid<end & N3_time_grid>start],start,end))
      for(i in 1:(length(grid_1_N3)-1)){
        temp<-data.frame(t=(grid_1_N3[i]+grid_1_N3[i+1])/2,tt1_less_t=first,tt1_less_t_delt1=second,tt1_less_t_tt1=third)
        temp2<-hal9001::make_design_matrix(as.matrix(temp),basis_list_N3_select)
        temp2<-c(1,as.vector(temp2))
        log_hazard_midpoint_N3<-sum(coef_N3*temp2)
        prob_N3_interval<-exp(-(grid_1_N3[i+1]-grid_1_N3[i])*exp(log_hazard_midpoint_N3))
        prob<-prob*prob_N3_interval
      }
    }
  }
  return(prob)
}


hazard_func<-function(t,first,second,third,basis_list_N3_select,coef_N3,type){
  if(type=="A1" |type=="N1"){
    temp<-data.frame(t=t,tt2_less_t=first,tt2_less_t_delt2=second,tt2_less_t_tt2=third)
  }else{
    temp<-data.frame(t=t,tt1_less_t=first,tt1_less_t_delt1=second,tt1_less_t_tt1=third)
  }

  temp2<-hal9001::make_design_matrix(as.matrix(temp),basis_list_N3_select)
  temp2<-c(1,as.vector(temp2))
  log_hazard_N3<-sum(coef_N3*temp2)
  hazard_N3<-exp(log_hazard_N3)
  return(hazard_N3)
}




