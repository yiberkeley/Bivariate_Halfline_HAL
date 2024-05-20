
#Specify the known conditional intensity when there is no censoring
#The conditional intensity specification for N1
set.seed(1000)
library(Rlab)
t1_gen<-runif(20,min=0,max=1)
t2_happened_gen<-rbern(20,0.5)
t2_happened_t2_gen<-t2_happened_gen*runif(20,min=0,max=1)
data_t1_basis_gen<-data.frame(t=t1_gen,t2_happened=t2_happened_gen,t2_happened_t2=t2_happened_t2_gen)
basis_list_t1 <- enumerate_basis(data_t1_basis_gen,max_degree = 2,smoothness_orders = 0)
basis_list_t1<-basis_list_t1[sample(1:length(basis_list_t1),floor(0.3*length(basis_list_t1)),replace=FALSE)]
coef_t1<-runif(length(basis_list_t1)+1,0,2)/4

#The conditional intensity specification for N2
t2_gen<-runif(20,min=0,max=1)
t1_happened_gen<-rbern(20,0.3)
t1_happened_t1_gen<-t1_happened_gen*runif(20,min=0,max=1)
data_t2_basis_gen<-data.frame(t=t2_gen,t1_happened=t1_happened_gen,t1_happened_t1=t1_happened_t1_gen)
basis_list_t2 <- enumerate_basis(data_t2_basis_gen,max_degree = 2,smoothness_orders = 0)
basis_list_t2<-basis_list_t2[sample(1:length(basis_list_t2),floor(0.3*length(basis_list_t2)),replace=FALSE)]
coef_t2<-runif(length(basis_list_t2)+1)/4

gen_t1_t2<-function(){
  grid_gen<-seq(0,1,by=0.001)

  t1_happened=0
  t1_happened_t1=0
  t2_happened=0
  t2_happened_t2=0
  for(i in 1:(length(grid_gen)-1)){

    if(t1_happened==0){
      #Jumping process N1
      temp<-data.frame(t=(grid_gen[i]+grid_gen[i+1])/2,t2_happened=t2_happened,t2_happened_t2=t2_happened_t2)
      temp2<-make_design_matrix(as.matrix(temp),basis_list_t1)
      temp2<-c(1,as.vector(temp2))
      log_hazard_midpoint_N1<-sum(coef_t1*temp2)
      prob_interval_N1<-exp(-(grid_gen[i+1]-grid_gen[i])*exp(log_hazard_midpoint_N1)) #Not jumping
      jump_interval_N1<-rbern(1,1-prob_interval_N1)
      if(jump_interval_N1==1){
        t1_happened=1
        t1_happened_t1=(grid_gen[i]+grid_gen[i+1])/2
      }
    }

    if(t2_happened==0){
      #Jumping process N2
      temp<-data.frame(t=(grid_gen[i]+grid_gen[i+1])/2,t1_happened=t1_happened ,t1_happened_t1=t1_happened_t1)
      temp2<-make_design_matrix(as.matrix(temp),basis_list_t2)
      temp2<-c(1,as.vector(temp2))
      log_hazard_midpoint_N2<-sum(coef_t2*temp2)
      prob_interval_N2<-exp(-(grid_gen[i+1]-grid_gen[i])*exp(log_hazard_midpoint_N2))
      jump_interval_N2<-rbern(1,1-prob_interval_N2)
      if(jump_interval_N2==1){
        t2_happened=1
        t2_happened_t2=(grid_gen[i]+grid_gen[i+1])/2
      }
    }
  }

  #if t1 hasn't happened, set to be 1
  if(t1_happened==0){
    t1=1
  }else{
    t1=t1_happened_t1
  }

  #if t2 hasn't happened, set to be 1
  if(t2_happened==0){
    t2=1
  }else{
    t2=t2_happened_t2
  }
  return(c(t1,t2))
}


datgen_knownIntensity_no_censoring<-function(n,seed){
  c1<-rep(4,n)
  c2<-rep(4,n)

  set.seed(seed)
  t1<-rep(0,n)
  t2<-rep(0,n)
  for(i in 1:n){
    t1_t2<-gen_t1_t2()
    t1[i]<-t1_t2[1]
    t2[i]<-t1_t2[2]
  }

  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)
  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}

datgen_knownIntensity_no_censoring_trunc<-function(n,seed){
  c1<-rep(4,n)
  c2<-rep(4,n)

  set.seed(seed)
  t1<-rep(0,n)
  t2<-rep(0,n)
  for(i in 1:n){
    t1_t2<-gen_t1_t2()
    t1[i]<-t1_t2[1]
    t2[i]<-t1_t2[2]
  }

  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)

  index_temp<-which(tt1>0.68)
  tt1[index_temp]<-0.68
  delt1[index_temp]<-0

  index_temp<-which(tt2>0.68)
  tt2[index_temp]<-0.68
  delt2[index_temp]<-0


  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}

