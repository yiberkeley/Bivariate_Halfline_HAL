
#flex_t helper function to get updated basis list and set up the penalty
flex_t_func<-function(data_obs,basis_list_N3,list_type,basis_option){
  #N3 jumping process
  #Get the poisitions of basis that only contain the t information only, change them to the grid
  #at failure times of corresponding process and set up penalty factor
  position_basis_N3<-lapply(basis_list_N3, `[[`, 1)
  position_basis_N3<-lapply(position_basis_N3,function(x){length(unlist(x))==1 & sum(unlist(x)==1)==1})
  position_basis_N3<-which(unlist(position_basis_N3)==TRUE)
  basis_list_N3_part_2<-basis_list_N3[-position_basis_N3]
  #Observed failure
  if(list_type=="N1"){
    N3_observed_failure_times<-(data_obs%>%filter(delt1==1))$tt1
  }else if(list_type=="A1"){
    N3_observed_failure_times<-(data_obs%>%filter(delt1==0))$tt1
  }else if(list_type=="N2"){
    N3_observed_failure_times<-(data_obs%>%filter(delt2==1))$tt2
  }else{
    N3_observed_failure_times<-(data_obs%>%filter(delt2==0))$tt2
  }
  N3_observed_failure_times<-N3_observed_failure_times[sample(1:length(N3_observed_failure_times),
                                                              ifelse(length(N3_observed_failure_times)<10,
                                                                     length(N3_observed_failure_times),
                                                                     floor(0.8*length(N3_observed_failure_times))), #previous 0.1
                                                              replace=FALSE)]
  if(basis_option=="0_order"){
    basis_list_N3_part_1 <- enumerate_basis(N3_observed_failure_times,max_degree = 2,smoothness_orders = 0)
  }else if(basis_option=="1_order"){
    basis_list_N3_part_1 <- enumerate_basis(N3_observed_failure_times,max_degree = 2,smoothness_orders = 1)
  }
  #Combined and set up penalty
  basis_list_N3<-c(basis_list_N3_part_1,basis_list_N3_part_2)
  penalty.factor_N3<-rep(1,length(basis_list_N3))
  penalty.factor_N3[1:length(basis_list_N3_part_1)]<-0.01
  result_list<-list()
  result_list[[1]]<-basis_list_N3
  result_list[[2]]<-penalty.factor_N3
  return(result_list)
}
