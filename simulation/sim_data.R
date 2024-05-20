library(data.table)

#Correlated high censoring
datgen_correlated_high_censoring<-function(n){
  temp<-rnorm2d(2*n, rho = 0.2) #old 0.7/0.2
  temp<-as.data.frame(temp)
  temp$V1<-temp$V1/6+0.5
  temp$V2<-temp$V2/6+0.5
  temp<-temp%>%filter(V1>0 & V2>0 &V1<1 &V2<1)
  sample_index<-sample(dim(temp)[1], n, replace = FALSE, prob = NULL)
  t1<-temp[sample_index,1]
  t2<-temp[sample_index,2]

  c1<-round(runif(n,0.0,1.2),4)
  c2<-round(runif(n,0.0,1.2),4)
  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)
  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}

#Correlated low censoring
datgen_correlated_low_censoring<-function(n){
  temp<-rnorm2d(2*n, rho = 0.2) #old 0.7
  temp<-as.data.frame(temp)
  temp$V1<-temp$V1/6+0.5
  temp$V2<-temp$V2/6+0.5
  temp<-temp%>%filter(V1>0 & V2>0 &V1<1 &V2<1)
  sample_index<-sample(dim(temp)[1], n, replace = FALSE, prob = NULL)
  t1<-temp[sample_index,1]
  t2<-temp[sample_index,2]

  c1<-round(runif(n,0.0,4),4)
  c2<-round(runif(n,0.0,4),4)
  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)
  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}

#Correlated no censoring
datgen_correlated_no_censoring<-function(n){
  temp<-rnorm2d(2*n, rho = 0.2) #old 0.7
  temp<-as.data.frame(temp)
  temp$V1<-temp$V1/6+0.5
  temp$V2<-temp$V2/6+0.5
  temp<-temp%>%filter(V1>0 & V2>0 &V1<1 &V2<1)
  sample_index<-sample(dim(temp)[1], n, replace = FALSE, prob = NULL)
  t1<-temp[sample_index,1]
  t2<-temp[sample_index,2]

  c1<-rep(4,n)
  c2<-rep(4,n)
  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)
  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}


#Independent low censoring
datgen_independent_low_censoring<-function(n){

  t1<-round(runif(n,0.0,1),4)
  t2<-round(runif(n,0.0,1),4)

  c1<-round(runif(n,0.0,4),4)
  c2<-round(runif(n,0.0,4),4)
  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)
  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}

#Independent high censoring
datgen_independent_high_censoring<-function(n){

  t1<-round(runif(n,0.0,1),4)
  t2<-round(runif(n,0.0,1),4)

  c1<-round(runif(n,0.0,1.2),4)
  c2<-round(runif(n,0.0,1.2),4)
  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)
  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}

#Independent no censoring
datgen_independent_no_censoring<-function(n){

  t1<-round(runif(n,0.0,1),4)
  t2<-round(runif(n,0.0,1),4)

  c1<-rep(4,n)
  c2<-rep(4,n)
  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)
  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}




#Correlated no censoring Trunc
datgen_correlated_no_censoring_trunc<-function(n){
  temp<-rnorm2d(2*n, rho = 0.2) #old 0.7
  temp<-as.data.frame(temp)
  temp$V1<-temp$V1/6+0.5
  temp$V2<-temp$V2/6+0.5
  temp<-temp%>%filter(V1>0 & V2>0 &V1<1 &V2<1)
  sample_index<-sample(dim(temp)[1], n, replace = FALSE, prob = NULL)
  t1<-temp[sample_index,1]
  t2<-temp[sample_index,2]

  c1<-rep(4,n)
  c2<-rep(4,n)
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



#Independent no censoring Trunc
datgen_independent_no_censoring_trunc<-function(n){

  t1<-round(runif(n,0.0,1),4)
  t2<-round(runif(n,0.0,1),4)

  c1<-rep(4,n)
  c2<-rep(4,n)
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

  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w)
  return(df)
}
