#flex_tail helper function to get updated basis list and set up the penalty
flex_tail_func<-function(data_obs,basis_list_N3,list_type,basis_option){
  #N3 jumping process
  #Get the positions of basis that contain t information and are larger than a threshold
  position_basis_N3_1<-lapply(basis_list_N3, `[[`, 1)
  position_basis_N3_1<-lapply(position_basis_N3_1,function(x){sum(unlist(x)==4)>=1 | sum(unlist(x)==1)>=1})
  position_basis_N3_1<-which(unlist(position_basis_N3_1)==TRUE)

  position_basis_N3_2<-lapply(basis_list_N3, `[[`, 2)
  position_basis_N3_2<-lapply(position_basis_N3_2,function(x){sum(unlist(x)>=0.7)>=1})
  position_basis_N3_2<-which(unlist(position_basis_N3_2)==TRUE)

  position_basis_N3<-intersect(position_basis_N3_1,position_basis_N3_2)

  #Change the penalty in those positions
  penalty.factor_N3<-rep(1,length(basis_list_N3))
  penalty.factor_N3[position_basis_N3_2]<-0.01

  result_list<-list()
  result_list[[1]]<-basis_list_N3
  result_list[[2]]<-penalty.factor_N3
  return(result_list)
}
