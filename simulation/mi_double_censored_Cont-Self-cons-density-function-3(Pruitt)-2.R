
# in this version i impute the half lines but do not remove the observed censored data point no half lines.
# the imputed ones do not contribute to p11.
# In p10 and p01, itr>=1 instead of itr==1. This means I am doing Pruitt's version where the prob of
# falling on a half line remains the same and do not get updated.


#rm(list=ls())

library(dplyr)
setwd("/Users/liyi/Documents/GitHub/Bivariate_Halfline_HAL/simulation")
source("utils.R")

datgen<-function(n){
  t1<-round(runif(n),4)
  t2<-round(runif(n),4)
  c1<-round(runif(n,0.0,2),4)
  c2<-round(runif(n,0.0,2),4)
  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)
  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}

norm_vec <- function(x){sqrt(crossprod(x))}

###START_CHANGE
draw10<-function(position,start,num_draws,res){
  temp_result<-Rejection_Sample_DensityBased(coef_N1_initial = res[[1]],
                                basis_list_N1_select = res[[2]],
                                coef_N2_initial = res[[3]],
                                basis_list_N2_select = res[[4]],
                                coef_A1_initial = res[[5]],
                                basis_list_A1_select = res[[6]],
                                coef_A2_initial =res[[7]],
                                basis_list_A2_select = res[[8]],
                                censoring=F,
                                half_line_position=position, #where the observed time is
                                half_line_type="N2", #which the censored event is
                                starting_of_half_line=start, # where the censored event happen
                                sample_size=num_draws
  )
  temp_result<-data.frame(values=temp_result[[1]],weights=temp_result[[2]],start=start,position=position)
  temp_result$weights<-temp_result$weights/sum(temp_result$weights) #Normalized weights
  temp_result$id<-seq(1:num_draws)
  return(temp_result)
}

draw01<-function(position,start,num_draws,res){
  temp_result<-Rejection_Sample_DensityBased(coef_N1_initial = res[[1]],
                                             basis_list_N1_select = res[[2]],
                                             coef_N2_initial = res[[3]],
                                             basis_list_N2_select = res[[4]],
                                             coef_A1_initial = res[[5]],
                                             basis_list_A1_select = res[[6]],
                                             coef_A2_initial =res[[7]],
                                             basis_list_A2_select = res[[8]],
                                             censoring=F,
                                             half_line_position=position, #where the observed time is
                                             half_line_type="N1", #which the censored event is
                                             starting_of_half_line=start, # where the censored event happen
                                             sample_size=num_draws
  )
  temp_result<-data.frame(values=temp_result[[1]],weights=temp_result[[2]],start=start,position=position)
  temp_result$weights<-temp_result$weights/sum(temp_result$weights) #Normalized weights
  temp_result$id<-seq(1:num_draws)
  return(temp_result)
}
###END_CHANGE

density_EM_bivar <- function(df,nrep){
  ###START_CHANGE
  #Run the hal fit
  res <- run_sim(n = n,
                 seed = seed,
                 df,
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
  ###END_CHANGE

  print("finished HAL fit")
  #set.seed(300);df<-datgen(50)

  df11<-filter(df, delt1==1&delt2==1)
  df10<-filter(df, delt1==1&delt2==0);dim(df10)
  df01<-filter(df, delt1==0&delt2==1)
  df00<-filter(df, delt1==0&delt2==0)

  df11$w<-1
  df00$w<-1

  #nrep<-5
  dup.dat10<-do.call("rbind", replicate(nrep, df10, simplify = FALSE))
  dup.dat01<-do.call("rbind", replicate(nrep, df01, simplify = FALSE))
  dup.dat00<-do.call("rbind", replicate(nrep, df00, simplify = FALSE))

  dup.rows<-duplicated(dup.dat10)
  ###START_CHANGE
  df_temp<-as.data.frame(cbind(df10$tt1,df10$tt2))
  draws10<-do.call(rbind,apply(df_temp,1,function(x){draw10(x[1],x[2],nrep-1,res)}))
  draws10<-draws10%>%arrange(id)
  dup.dat10$tt2[dup.rows]<-draws10$values
  dup.dat10$delt2[dup.rows]<-1
  dup.dat10$type[dup.rows]<-0
  dup.dat10$w[dup.rows]<-draws10$weights
  ###END_CHANGE

  dup.rows<-duplicated(dup.dat01)
  ###START_CHANGE
  df_temp<-as.data.frame(cbind(df01$tt2,df01$tt1))
  draws01<-do.call(rbind,apply(df_temp,1,function(x){draw01(x[1],x[2],nrep-1,res)}))
  draws01<-draws01%>%arrange(id)
  dup.dat01$tt1[dup.rows]<-draws01$values
  dup.dat01$delt1[dup.rows]<-1
  dup.dat01$type[dup.rows]<-0
  dup.dat01$w[dup.rows]<-draws01$weights
  ###END_CHANGE

  dup.rows<-duplicated(dup.dat00)
  dup.dat00$tt1[dup.rows]<-runif(sum(dup.rows),dup.dat00$tt1[dup.rows] ,1)
  dup.dat00$delt1[dup.rows]<-1
  dup.dat00$tt2[dup.rows]<-runif(sum(dup.rows),dup.dat00$tt2[dup.rows] ,1)
  dup.dat00$delt2[dup.rows]<-1
  dup.dat00$type[dup.rows]<-0
  dup.dat00$w[dup.rows]<-1/(nrep-1)

  print("finished augmentation")

  aug.df<-rbind(df11,dup.dat10,dup.dat01,dup.dat00);
  full.df<-aug.df%>%filter(delt1==1&delt2==1)
  aug.df<-rbind(df11,dup.dat10,dup.dat01,df00)

  observed.f.func<-matrix(NA,ncol=10,nrow=nrow(aug.df))
  observed.true.s.func<-observed.s.func<-rep(0,nrow(aug.df))

  #  plot1<-ggplot(df,aes(x=tt1, y=tt2, shape=as.factor(delt2),color=as.factor(delt2))) +
  #    geom_point(alpha = 0.7, size = 2) +scale_shape_manual(values = c(17,19))+ labs(shape = "delta 2",color= "delta 2")

  #  plot1<-ggplot(temp.df10,aes(x=tt1, y=tt2, shape=as.factor(delt2),color=as.factor(delt2))) +
  #    geom_point(alpha = 0.7, size = 2) +scale_shape_manual(values = c(17,19))+ labs(shape = "delta 2",color= "delta 2")

  big.grid.t<-round((cbind(rep(seq(0.0,1,0.05),each=21),rep(seq(0,1,0.05),21))),2)
  #big.grid.t<-round(select(aug.df,c("tt1","tt2")),2)
  s.true<-(1-big.grid.t[,1])*(1-big.grid.t[,2])
  f.func<-s.func<-matrix(NA,nrow=nrow(big.grid.t),ncol=2)


  for (itr in 1:5){

    for(ii in 1:nrow(aug.df)){
      #      for(ii in 1:10){

      grid.t<-cbind(aug.df[ii,c("tt1")],aug.df[ii,c("tt2")])
      delt<-cbind(aug.df[ii,c("delt1")],aug.df[ii,c("delt2")])
      w<-aug.df[ii,"w"]

      p11<-sum((aug.df$type==1)&(aug.df$delt1==1)&(aug.df$delt2==1)&(aug.df$tt1==grid.t[1])&(aug.df$tt2==grid.t[2]))


      temp.df10<- df %>% filter(delt1==1,delt2==0,tt1==grid.t[1])
      temp.df01<- df %>% filter(delt1==0,delt2==1,tt2==grid.t[2])
      temp.df00<- df %>% filter(delt1==0,delt2==0,tt2<=grid.t[2],tt1<=grid.t[1])

      p10<-0
      if(nrow(temp.df10)!=0) {
        for(iii in 1:nrow(temp.df10)){
          co<-temp.df10[iii,c("tt1","tt2")]
          crit<-(aug.df$tt1==grid.t[1])&(aug.df$tt2==grid.t[2])
          if(itr>=1){
            cond.prob10<-(as.numeric(co[2])<grid.t[2])*1/
              sum((aug.df$tt1==grid.t[1])&(aug.df$tt2>as.numeric(co[2])))
          }else{
            cond.prob10<-(as.numeric(co[2])<grid.t[2])*sum ((observed.f.func[crit,itr-1])/
                                                              sum((observed.f.func[(aug.df$tt1==grid.t[1])&(aug.df$tt2>as.numeric(co[2])),itr-1]) ) )}
          temp.w10<- cond.prob10
          p10<-temp.w10+p10
        }
      }

      p01<-0
      if(nrow(temp.df01)!=0) {
        for(iii in 1:nrow(temp.df01)){
          co<-temp.df01[iii,c("tt1","tt2")]
          crit<-(aug.df$tt1==grid.t[1])&(aug.df$tt2==grid.t[2])
          if(itr>=1){
            cond.prob01<-(as.numeric(co[1])<grid.t[1])*1/
              sum((aug.df$tt2==grid.t[2])&(aug.df$tt1>as.numeric(co[1])))

          }else{
            cond.prob01<-(as.numeric(co[1])<grid.t[1])*sum ((observed.f.func[crit,itr-1])/
                                                              sum((observed.f.func[(aug.df$tt2==grid.t[2])&(aug.df$tt1>as.numeric(co[1])),itr-1]) ) )}
          temp.w01<- cond.prob01
          p01<-temp.w01+p01
        }
      }


      p00<-0
      if(nrow(temp.df00)!=0) {
        for(iii in 1:nrow(temp.df00)){
          co<-temp.df00[iii,c("tt1","tt2")]
          crit<-(aug.df$tt1==grid.t[1])&(aug.df$tt2==grid.t[2])
          if(itr==1){
            if(nrep>2){  cond.prob00<-(as.numeric(co[1])<grid.t[1])*(as.numeric(co[2])<grid.t[2])*delt[1]*delt[2]*w/
              sum(aug.df$w[(aug.df$tt2>as.numeric(co[2]))&(aug.df$tt1>as.numeric(co[1]))]) }else{
                cond.prob00<-(as.numeric(co[1])<grid.t[1])*(as.numeric(co[2])<grid.t[2])*delt[1]*delt[2]*1/
                  sum((aug.df$tt2>as.numeric(co[2]))&(aug.df$tt1>as.numeric(co[1]))) }

          }else{
            cond.prob00<-(as.numeric(co[1])<grid.t[1])*(as.numeric(co[2])<grid.t[2])*delt[1]*delt[2]*
              sum ((observed.f.func[crit,itr-1])/
                     sum((observed.f.func[(aug.df$tt2>as.numeric(co[2]))&(aug.df$tt1>as.numeric(co[1])),itr-1]) ) )}
          temp.w00<- cond.prob00;#print(cond.prob00)
          p00<-temp.w00+p00
        }

      }

  #    print(c(ii,p11,p10,p01,p00))


      #   if(itr==1){observed.f.func[ii,1]<-runif(1)}else{
      observed.f.func[ii,itr]<- (p11+p01+p10+1 *p00)*1/nrow(df)
      #                                                  }

    }
  }

  #observed.f.func[,5]
  #  s.dat<-as.data.frame((cbind(aug.df[,c("tt1","tt2")]
  #                                   ,observed.f.func[,1],observed.f.func[,4],1/11^2)) )

    observed.f.func[,c(1,2,5)]

  # if i geenrate more than one replicate then f converges to different values
  # depending on the initial but the s functions are identical.

  ##  for(i in 1:nrow(aug.df)){
  ##    grid.t<-cbind(aug.df[i,c("tt1")],aug.df[i,c("tt2")])
  ##    row<-which(grid.t[1]<aug.df$tt1 & grid.t[2]<aug.df$tt2)
  ##    observed.s.func[i]<-sum(observed.f.func[row,itr])
  ##    observed.true.s.func[i]<-(1-grid.t[1])*(1-grid.t[2])
  ##  }

  ##  cbind(observed.true.s.func,observed.s.func)


  for(i in 1:nrow(big.grid.t)){
    grid.t<-big.grid.t[i,]

    #    row<-which(big.grid.t[,1]==round(aug.df$tt1[i],2) & big.grid.t[,2]==round(aug.df$tt2[i],2))
    #    f.func[row,1]<-observed.f.func[i,5]/observed.f.func[i,1]
    row<-which(big.grid.t[i,1]<aug.df$tt1 & big.grid.t[i,2]<aug.df$tt2)
    s.func[i,1]<-sum(observed.f.func[row,itr])
    s.func[i,2]<-sum(observed.f.func[row,1])
  }
    s.true<-(1-big.grid.t[,1])*(1-big.grid.t[,2])
  ##      norm_vec(s.func[,1]-s.true)
  ##      norm_vec(s.func[,2]-s.true)
  ##  cbind(s.func[,1],s.true)

  ##########################################################################
  ##########################################################################
  ##########################################################################
  ##########################################################################
#######  s.true<-NULL
#  full.df<-df
#  df11<-filter(df, delt1==1&delt2==1)
#  full.df10<-filter(df, delt1==1&delt2==0);dim(df10)
  # full.df01<-filter(df, delt1==0&delt2==1)
  # full.df00<-filter(df, delt1==0&delt2==0)
  #
  #
  # full.df10$tt2<-runif(nrow(df10),df10$tt2 ,1)
  # full.df01$tt1<-runif(nrow(df01),df01$tt1 ,1)
  # full.df00$tt1<-runif(nrow(df00),df00$tt1 ,1)
  # full.df00$tt2<-runif(nrow(df00),df00$tt2 ,1)
  # full.df<-rbind(df11,full.df10,full.df01,full.df00)

  for(i in 1:nrow(big.grid.t)){
    grid.t<-big.grid.t[i,]
    row<-(big.grid.t[i,1]<full.df$tt1 & big.grid.t[i,2]<full.df$tt2)
    ########s.true[i]<-sum(row)/nrow(full.df)
  }
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ##########################################################################


  #   s.dat%>% filter(tt1==.7 & tt2==.4)

  #result<-observed.f.func[,1]#round(cbind(df[,c("tt1","tt2")] ,observed.f.func[,1],observed.f.func[,5],1/11^2),3)
#  result<-s.func[,1]#norm_vec(s.func[,1]-s.true)*sqrt(nrow(df))
#  result<-list(estimated.survival=s.func[,1],true.survival=s.true)
  result<-c(estimated.survival=s.func[,1],true.survival=s.true)
  return(result)
}

big.grid.t<-round((cbind(rep(seq(0.0,1,0.05),each=21),rep(seq(0,1,0.05),21))),2)
s.true<-(1-big.grid.t[,1])*(1-big.grid.t[,2])
norm_vec <- function(x){sqrt(crossprod(x))}


n.vec<-c(seq(520,980,by=20))
n.vec<-100#seq(10,1950,by=50)
nrep=4
scaled.l2norm<-scaled.bias<-NULL
for(ii in 1:length(n.vec)){

  n<-n.vec[ii]
  nsim<-10#n/10
sim.1 <- vector("list", length = nsim)

set.seed(320)
for(i in 1:nsim){
  dat<-datgen(n = n )
  sim.1[[i]]<-dat
}

library(doParallel)#Load parallel library
no_cores <- detectCores()-1 #Set number of cores
cl <- makeCluster(no_cores)
registerDoParallel((cl))
###START_CHANGE
# results.list <- foreach(k =1:length(sim.1), .packages = c("MASS","dplyr")) %dopar%
#   tryCatch(density_EM_bivar(df = sim.1[[k]],nrep), error = function(e) {print(e); NA})
results.list<-lapply(1:length(sim.1),function(x){density_EM_bivar(df = sim.1[[x]],nrep)})
###END_CHANGE
stopCluster(cl)
temp<-as.data.frame(do.call(rbind,results.list))
true.survival.df<-temp%>%dplyr::select(starts_with("true"))
est.survival.df<-temp%>%dplyr::select(starts_with("est"))

true.survival<-apply(true.survival.df,2,function(x) mean(x,na.rm=TRUE))
#true.survival
est.survival<-apply(est.survival.df,2,function(x) mean(x,na.rm=TRUE))
#est.survival
bias.survival<-apply(est.survival.df-true.survival.df,2,function(x) mean(x,na.rm=TRUE))

scaled.bias[ii]<-sum(abs((est.survival-true.survival)*sqrt(n)))/nrow(big.grid.t)

norm_vec <- function(x){sqrt(crossprod(x))}
scaled.l2norm[ii]<-norm_vec(est.survival-true.survival)*sqrt(n)

note<-paste("I generated nrep points on each half line and double censored.
            i approximated true using full data.
            the distribution over half liness does not change just like
            Pruitt as it does not change the density on half lines.
            nsim=0.1n")
 save(note,nrep,true.survival.df,est.survival.df,est.survival,true.survival,
      results.list,big.grid.t,n,nsim,density_EM_bivar,datgen,
      file = paste0("BivSurv_IndepUnif_sim_V3_Pruitt_n=",n,"nrep=",nrep,"nsim=",nsim,".rda"))

print(n)
}

note<-paste(" in this version i impute the half lines but do not remove the observed censored data point no half lines.
 the imputed ones do not contribute to p11.
 In p10 and p01, itr>=1 instead of itr==1. This means I am doing Pruitt's version where the prob of
 falling on a half line remains the same and do not get updated.
            nsim=0.1n")
#save(note,nrep,scaled.bias,scaled.l2norm,big.grid.t,n.vec,density_EM_bivar,datgen,
#     file = paste0("BivSurv_IndepUnif_all_nvec_sim_n=",n,"nrep=",nrep,"nsim=",nsim,".rda"))




