
library(ggplot2)
library(ggpubr)
library(purrr)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)
library(future.apply)
library(fMultivar)
source("sim_data.R")

`%+%` <- function(a, b) paste0(a, b)

#when n=500,  perc_sample=0.1, perc_sample_basis=0.1
#when n=1000,  perc_sample=0.02, perc_sample_basis=0.02
#when n=1000,  perc_sample=0.01, perc_sample_basis=0.01

#In the lambda_flex, nearly half basis selected are the failure time based indicators to modeling the t part.
#In the lambda case, only 6 percent selected are the failure time based indicators.

#Evaluate the truth
datgen<-function(n){return(datgen_correlated_no_censoring(n))} #change here
t1_plot<- seq(0,1,0.1)
t2_plot <- seq(0,1,0.1)
# t1_plot<- t1.detail
# t2_plot <- t2.detail
plot_data<-expand.grid(t1=t1_plot,t2=t2_plot)
set.seed(100)
nCores<-detectCores()-2
plan(multisession,workers = nCores)
data_full_sim<-datgen(10000)
truth<-future_apply(plot_data,1,function(x) mean(data_full_sim$t1>x[1] & data_full_sim$t2>x[2]))

truth.0.2corr.nc<-truth
truth.0.2corr.lc<-truth
truth.0.2corr.hc<-truth
truth.indep.nc<-truth
truth.indep.lc<-truth
truth.indep.hc<-truth

truth.0.7corr.nc<-truth
truth.0.7corr.lc<-truth
truth.0.7corr.hc<-truth

#Get the empirical naive
set.seed(101)
data_full_sample<-datgen(500)
data_obs_sample<-data_full_sample%>%select(tt1,delt1,tt2,delt2)
empirical_truth<-future_apply(plot_data,1,function(x) mean(data_obs_sample$tt1>x[1] & data_obs_sample$tt2>x[2]))

#Input the simulation result
input_result_fun<-function(filename_first,seed_vals){
  result<-list()
  for(i in 1:length(seed_vals)){
    result_temp<-readRDS(filename_first %+% seed_vals[i] %+% ".RDS")
    result[[i]]<-result_temp
  }
  return(result)
}

# result.l1<-input_result_fun(filename_first="cl_022924/cl_022624_usual_0_order_masking_rightpoint4_500_",
#                            seed_vals = c(101,500,2000,4321))

# result.l1<-input_result_fun(filename_first="cl_022224/cl_022324_usual_rightpoint4_500_",
#                             seed_vals = c(101,500,2000,4321,10113,48001))

#For high correlation 0.7, low censoring, true propensity
result.l1.1<-input_result_fun(filename_first="cl_022924/cl_030424_usual_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                            seed_vals = c(101,500,2000,4321))


result.l1.2<-input_result_fun(filename_first="cl_022924/cl_030424_flex_t_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001))

result.l2.1<-input_result_fun(filename_first="cl_022924/cl_030424_usual_0_order_tailShrinkage_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321,8688))

result.l2.2<-input_result_fun(filename_first="cl_022924/cl_030424_flex_t_0_order_tailShrinkage_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001,50002))

result.l3.1<-input_result_fun(filename_first="cl_022924/cl_030424_usual_0_order_masking_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321))

result.l3.2<-input_result_fun(filename_first="cl_022924/cl_022924_flex_t_0_order_masking_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001,50002))




#For low correlation 0.2, low censoring, true propensity
result.l1.1.low<-input_result_fun(filename_first="cl_030524_0.2corr/cl_030524_0.2corr_usual_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321,8688,10113))


result.l1.2.low<-input_result_fun(filename_first="cl_030524_0.2corr/cl_030524_0.2corr_flex_t_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000))

result.l2.1.low<-input_result_fun(filename_first="cl_030524_0.2corr/cl_030524_0.2corr_usual_0_order_tailShrinkage_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001))

result.l2.2.low<-input_result_fun(filename_first="cl_030524_0.2corr/cl_030524_0.2corr_flex_t_0_order_tailShrinkage_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001))

result.l3.1.low<-input_result_fun(filename_first="cl_030524_0.2corr/cl_030524_0.2corr_usual_0_order_masking_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001))

result.l3.2.low<-input_result_fun(filename_first="cl_030524_0.2corr/cl_030524_0.2corr_flex_t_0_order_masking_rightpoint4_rightSubsampleAndPenalty_500_",
                              seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001,50002))


#For independent , low censoring, true propensity
result.indep.low.1<-input_result_fun(filename_first="030624_various/cl_030624_independent_usual_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                                  seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001))

#For independent , high censoring, true propensity
result.indep.high.1<-input_result_fun(filename_first="030624_various/ch_030624_independent_usual_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                                   seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001))

#For low correlation 0.2, high censoring, true propensity
result.corr0.2.high.1<-input_result_fun(filename_first="030624_various/ch_030624_0.2corr_usual_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                                   seed_vals = c(101,500,2000))

#For high correlation 0.7, high censoring, true propensity
result.corr0.7.high.1<-input_result_fun(filename_first="030624_various/ch_030624_0.7corr_usual_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                                   seed_vals = c(101,500,2000,4321,8688,10113))


#For low correlation, low censoring, true propensity with undersmoothing
result.l4.1.low<-input_result_fun(filename_first="030824_various/cl_030824_0.2corr_usual_undersmooth_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                                  seed_vals = c(101,500,2000,4321,8688,10113))
result.l4.2.low<-input_result_fun(filename_first="030824_various/cl_030824_0.2corr_flex_t_undersmooth_0_order_rightpoint4_rightSubsampleAndPenalty_500_",
                                  seed_vals = c(101,500,2000,4321,8688,10113))

##########
#Increase the sample size, test the undersmoothing.
#########

#####
#One step estimator (-PnD^*)
#####

#In full, 1-\hat{P} solved or not, greater than the 1/\sqrt(n)
#Does the TMLE equal to the empirical.




#For full data
result.full.indep<-input_result_fun(filename_first="030824_various/cn_031224_independent_usual_full_0_order_rightpoint4_rightSubsampleAndPenalty_100_",
                                  seed_vals = c(101,500,2000,4321,8688,10113))
result.full.corr0.2<-input_result_fun(filename_first="030824_various/cn_031224_0.2corr_usual_full_0_order_rightpoint4_rightSubsampleAndPenalty_100_",
                                  seed_vals = c(101,500,2000,4321,8688,10113))

#For full data (more detailed version)
result.full.indep.detail<-input_result_fun(filename_first="030824_various/cn_031324_independent_usual_full_0_order_rightpoint4_rightSubsampleAndPenalty_rightIntegration100_",
                                    seed_vals = c(101,500,2000,4321,8688,10113,36899,45000))
result.full.corr0.2.detail<-input_result_fun(filename_first="030824_various/cn_031324_0.2corr_usual_full_0_order_rightpoint4_rightSubsampleAndPenalty_rightIntegration100_",
                                      seed_vals = c(101,500,2000,4321,8688,10113,36899,45000))


t1.detail<-result.full.indep.detail[[1]][[1]][[1]]$t1
t2.detail<-result.full.indep.detail[[1]][[1]][[1]]$t2


t1<-result.l1.1[[1]][[1]][[1]]$t1
t2<-result.l1.1[[1]][[1]][[1]]$t2


##Evaluate over the whole grid
#Create bias_sd_grid_data
bias_sd_grid_fun<-function(result_list,truth){
  bias_sd<-NULL
  for(j in 1:length(result_list)){
      bias_temp<-mean(abs(pmax(result_list[[j]][[1]][[1]]$survival_curve-result_list[[j]][[1]][[1]]$survival_curve[1]+1,0)-truth))
      sd_temp<-sd(pmax(result_list[[j]][[1]][[1]]$survival_curve-result_list[[j]][[1]][[1]]$survival_curve[1]+1,0)-truth)
      N1_l1_temp<-result_list[[j]][[1]][[2]]
      N2_l1_temp<-result_list[[j]][[1]][[3]]
      bias_sd<-rbind(bias_sd,c(N1_l1_temp,N2_l1_temp,bias_temp,sd_temp,update=F))
  }

  for(j in 1:length(result_list)){
    bias_temp<-mean(abs(pmax(result_list[[j]][[2]][[1]]$survival_curve-result_list[[j]][[2]][[1]]$survival_curve[1]+1,0)-truth))
    sd_temp<-sd(pmax(result_list[[j]][[2]][[1]]$survival_curve-result_list[[j]][[2]][[1]]$survival_curve[1]+1,0)-truth)
    N1_l1_temp<-result_list[[j]][[2]][[2]]
    N2_l1_temp<-result_list[[j]][[2]][[3]]
    bias_sd<-rbind(bias_sd,c(N1_l1_temp,N2_l1_temp,bias_temp,sd_temp,update=T))
  }

  return(bias_sd)
}

result_temp<-result.l1.8
bias_sd_grid_data_l1<-bias_sd_grid_fun(result_temp,truth = truth)
colnames(bias_sd_grid_data_l1)<-c("N1_l1","N2_l1","bias","sd","update")
bias_sd_grid_data_l1<-as.data.frame(bias_sd_grid_data_l1)

ggplot(bias_sd_grid_data_l1,aes(x=as.factor(update),y=bias))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  ggtitle("avg|residual|")

ggplot(bias_sd_grid_data_l1,aes(x=as.factor(update),y=sd))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  ggtitle("sd(residual)")


##################
##################
##################
##################

#Investigation into the fit

two_plot_data_func<-function(result1,result2,t1,t2,truth1,truth2){
  result_temp<-result1
  sum_surve_initial<-  pmax(result_temp[[1]][[1]][[1]]$survival_curve-result_temp[[1]][[1]][[1]]$survival_curve[1]+1,0)
  sum_surve_update<- pmax(result_temp[[1]][[2]][[1]]$survival_curve-result_temp[[1]][[2]][[1]]$survival_curve[1]+1,0)
  for(i in 2: length(result_temp)){
    sum_surve_initial<-sum_surve_initial+pmax(result_temp[[i]][[1]][[1]]$survival_curve-result_temp[[i]][[1]][[1]]$survival_curve[1]+1,0)
    sum_surve_update<-sum_surve_update+pmax(result_temp[[i]][[2]][[1]]$survival_curve-result_temp[[i]][[2]][[1]]$survival_curve[1]+1,0)
  }
  mean_surve_initial1<-sum_surve_initial/length(result_temp)
  mean_surve_update1<-sum_surve_update/length(result_temp)

  rss_surve_initial<-(pmax(result_temp[[1]][[1]][[1]]$survival_curve-result_temp[[1]][[1]][[1]]$survival_curve[1]+1,0)-mean_surve_initial1)^2
  rss_surve_update<-(pmax(result_temp[[1]][[2]][[1]]$survival_curve-result_temp[[1]][[2]][[1]]$survival_curve[1]+1,0)-mean_surve_update1)^2
  for(i in 2: length(result_temp)){
    rss_surve_initial<-rss_surve_initial+(pmax(result_temp[[i]][[1]][[1]]$survival_curve-result_temp[[i]][[1]][[1]]$survival_curve[1]+1,0)-mean_surve_initial1)^2
    rss_surve_update<-rss_surve_update+(pmax(result_temp[[i]][[2]][[1]]$survival_curve-result_temp[[i]][[2]][[1]]$survival_curve[1]+1,0)-mean_surve_update1)^2
  }
  sd_surve_initial1<-sqrt(rss_surve_initial/length(result_temp))
  sd_surve_update1<-sqrt(rss_surve_update/length(result_temp))



  plot_data1.1<-data_frame(t1=t1,t2=t2,mean_surve_bias=mean_surve_initial1-truth1, sd_surve=sd_surve_initial1, type=1,update=0)
  plot_data1.2<-data_frame(t1=t1,t2=t2,mean_surve_bias=mean_surve_update1-truth2 , sd_surve=sd_surve_update1,  type=1,update=1)

  result_temp<-result2
  sum_surve_initial<- pmax(result_temp[[1]][[1]][[1]]$survival_curve-result_temp[[1]][[1]][[1]]$survival_curve[1]+1,0)
  sum_surve_update<-pmax(result_temp[[1]][[2]][[1]]$survival_curve-result_temp[[1]][[2]][[1]]$survival_curve[1]+1,0)
  for(i in 2: length(result_temp)){
    sum_surve_initial<-sum_surve_initial+ pmax(result_temp[[i]][[1]][[1]]$survival_curve-result_temp[[i]][[1]][[1]]$survival_curve[1]+1,0)
    sum_surve_update<-sum_surve_update+pmax(result_temp[[i]][[2]][[1]]$survival_curve-result_temp[[i]][[2]][[1]]$survival_curve[1]+1,0)
  }
  mean_surve_initial2<-sum_surve_initial/length(result_temp)
  mean_surve_update2<-sum_surve_update/length(result_temp)

  rss_surve_initial<-(pmax(result_temp[[1]][[1]][[1]]$survival_curve-result_temp[[1]][[1]][[1]]$survival_curve[1]+1,0)-mean_surve_initial2)^2
  rss_surve_update<-(pmax(result_temp[[1]][[2]][[1]]$survival_curve-result_temp[[1]][[2]][[1]]$survival_curve[1]+1,0)-mean_surve_update2)^2
  for(i in 2: length(result_temp)){
    rss_surve_initial<-rss_surve_initial+(pmax(result_temp[[i]][[1]][[1]]$survival_curve-result_temp[[i]][[1]][[1]]$survival_curve[1]+1,0)-mean_surve_initial1)^2
    rss_surve_update<-rss_surve_update+(pmax(result_temp[[i]][[2]][[1]]$survival_curve-result_temp[[i]][[2]][[1]]$survival_curve[1]+1,0)-mean_surve_update1)^2
  }
  sd_surve_initial2<-sqrt(rss_surve_initial/length(result_temp))
  sd_surve_update2<-sqrt(rss_surve_update/length(result_temp))



  plot_data2.1<-data_frame(t1=t1,t2=t2,mean_surve_bias=mean_surve_initial2-truth1, sd_surve=sd_surve_initial2, type=2,update=0)
  plot_data2.2<-data_frame(t1=t1,t2=t2,mean_surve_bias=mean_surve_update2-truth2 , sd_surve=sd_surve_update2,  type=2,update=1)

  return(rbind(plot_data1.1,plot_data1.2,plot_data2.1,plot_data2.2))
}

#Any two results
result1<-result.indep.low.1 #independent low censoring
truth1<-truth.indep.lc
result2<-result.l1.1.low #0.2 corr low censoring
truth2<-truth.0.2corr.lc


result1<-result.indep.high.1 #independent high censoring
truth1<-truth.indep.hc
result2<-result.corr0.2.high.1 #0.2 corr high censoring
truth2<-truth.0.2corr.hc

result1<-result.l4.1.low #0.2 corr, low censoring, undersmooth,  usual
truth1<-truth.0.2corr.lc
result2<-result.l4.2.low #0.2 corr, low censoring, undersmooth,  flex t
truth2<-truth.0.2corr.lc


#

#########
#########
#CEHCK THE 1()-s().., Also check the integration one to see whether there is error.
#TMLE equals survival function. four point should be empirical survival fucntion at these four points. (one point)
#########
#########
#MAKE THE SUM OF THE TAIL OF THE DENSITY EQUALS TO THAT OF THE EMPIRICAL

result1<-result.full.indep #independent, full
truth1<-truth.indep.nc
result2<-result.full.corr0.2 #0.2 corr, full
truth2<-truth.0.2corr.nc

# result1<-result.l1.1 #0.7 corr, low censroing
# truth1<-truth.0.7corr.lc
# result2<-result.corr0.7.high.1 #0.7 corr, high censoring
# truth2<-truth.0.7corr.hc


library(plotly)
plot_data<-two_plot_data_func(result1,result2,t1,t2,truth1,truth2)
plot_ly(plot_data%>%filter(update==0),x=~t1,y=~t2,z=~mean_surve_bias,color=~type,colors = c('#BF382A', '#0C4B8E')) #initial
plot_ly(plot_data%>%filter(update==1),x=~t1,y=~t2,z=~mean_surve_bias,color=~type,colors = c('#BF382A', '#0C4B8E')) #update

plot_ly(plot_data%>%filter(type==1),x=~t1,y=~t2,z=~mean_surve_bias,color=~update,colors = c("#009999", "#0000FF")) #result1
plot_ly(plot_data%>%filter(type==2),x=~t1,y=~t2,z=~mean_surve_bias,color=~update,colors = c("#009999", "#0000FF")) #result2


#Look into good region
plot_ly(plot_data%>%filter(update==0 & t1<=0.7 & t2<=0.7),x=~t1,y=~t2,z=~mean_surve_bias,color=~type,colors = c('#BF382A', '#0C4B8E')) #initial
plot_ly(plot_data%>%filter(update==1 & t1<=0.7 & t2<=0.7),x=~t1,y=~t2,z=~mean_surve_bias,color=~type,colors = c('#BF382A', '#0C4B8E')) #update

plot_ly(plot_data%>%filter(type==1 & t1<=0.7 & t2<=0.7),x=~t1,y=~t2,z=~mean_surve_bias,color=~update,colors = c("#009999", "#0000FF")) #result1
plot_ly(plot_data%>%filter(type==2 & t1<=0.7 & t2<=0.7),x=~t1,y=~t2,z=~mean_surve_bias,color=~update,colors = c("#009999", "#0000FF")) #result2


#One specific fit
type_look <-1
update_look<-1
temp<-plot_data%>%filter(update==update_look,type==type_look)%>%select(-update,-type)
long <- melt(setDT(temp), id.vars = c("t1","t2"), variable.name = "type")
plot_ly(long,x=~t1,y=~t2,z=~value,color=~type,colors = c("magenta2", "orchid4")) #bias vs sd

long2 <- temp%>%mutate(bias_sd_ratio=ifelse(sd_surve!=0,mean_surve_bias/sd_surve,0))
plot_ly(long2,x=~t1,y=~t2,z=~abs(bias_sd_ratio)) # bias/sd
hist(abs(long2$bias_sd_ratio))

summary(abs(long2$bias_sd_ratio))


##################
##################
##################
##################

#Merge two fit

merge_plot_data_func<-function(result1,result2,t1,t2,truth){
  result_temp<-result1
  result_temp2<-result2
  sum_surve_initial<-(result_temp[[1]][[1]][[1]]$survival_curve/result_temp[[1]][[1]][[1]]$survival_curve[1]+
    result_temp2[[1]][[1]][[1]]$survival_curve/result_temp2[[1]][[1]][[1]]$survival_curve[1])/2

  sum_surve_update<-(result_temp[[1]][[2]][[1]]$survival_curve/result_temp[[1]][[2]][[1]]$survival_curve[1]+
    result_temp2[[1]][[2]][[1]]$survival_curve/result_temp2[[1]][[2]][[1]]$survival_curve[1])/2

  for(i in 2: length(result_temp)){
    sum_surve_initial<-sum_surve_initial+(result_temp[[i]][[1]][[1]]$survival_curve/result_temp[[i]][[1]][[1]]$survival_curve[1]+
                                            result_temp2[[i]][[1]][[1]]$survival_curve/result_temp2[[i]][[1]][[1]]$survival_curve[1])/2

    sum_surve_update<-sum_surve_update+(result_temp[[i]][[2]][[1]]$survival_curve/result_temp[[i]][[2]][[1]]$survival_curve[1]+
                                          result_temp2[[i]][[2]][[1]]$survival_curve/result_temp2[[i]][[2]][[1]]$survival_curve[1])/2
  }
  mean_surve_initial1<-sum_surve_initial/length(result_temp)
  mean_surve_update1<-sum_surve_update/length(result_temp)

  rss_surve_initial<-((result_temp[[1]][[1]][[1]]$survival_curve/result_temp[[1]][[1]][[1]]$survival_curve[1]+
                         result_temp2[[1]][[1]][[1]]$survival_curve/result_temp2[[1]][[1]][[1]]$survival_curve[1])/2-mean_surve_initial1)^2

  rss_surve_update<-((result_temp[[1]][[2]][[1]]$survival_curve/result_temp[[1]][[2]][[1]]$survival_curve[1]+
                        result_temp2[[1]][[2]][[1]]$survival_curve/result_temp2[[1]][[2]][[1]]$survival_curve[1])/2-mean_surve_update1)^2

  for(i in 2: length(result_temp)){
    rss_surve_initial<-rss_surve_initial+((result_temp[[i]][[1]][[1]]$survival_curve/result_temp[[i]][[1]][[1]]$survival_curve[1]+
                                             result_temp2[[i]][[1]][[1]]$survival_curve/result_temp2[[i]][[1]][[1]]$survival_curve[1])/2-mean_surve_initial1)^2

    rss_surve_update<-rss_surve_update+((result_temp[[i]][[2]][[1]]$survival_curve/result_temp[[i]][[2]][[1]]$survival_curve[1]+
                                           result_temp2[[i]][[2]][[1]]$survival_curve/result_temp2[[i]][[2]][[1]]$survival_curve[1])/2-mean_surve_update1)^2
  }
  sd_surve_initial1<-sqrt(rss_surve_initial/length(result_temp))
  sd_surve_update1<-sqrt(rss_surve_update/length(result_temp))



  plot_data1.1<-data_frame(t1=t1,t2=t2,mean_surve_bias=mean_surve_initial1-truth, sd_surve=sd_surve_initial1,update=0)
  plot_data1.2<-data_frame(t1=t1,t2=t2,mean_surve_bias=mean_surve_update1-truth , sd_surve=sd_surve_update1 ,update=1)

  return(rbind(plot_data1.1,plot_data1.2))
}

#Look into masking update and regular update
result1<-result.l1.1
result2<-result.l3.1

library(plotly)
plot_data<-merge_plot_data_func(result1,result2,t1,t2,truth)

update_look<-0
temp<-plot_data%>%filter(update==update_look)%>%select(-update)
long <- melt(setDT(temp), id.vars = c("t1","t2"), variable.name = "type")
plot_ly(long,x=~t1,y=~t2,z=~value,color=~type,colors = c("magenta2", "orchid4")) #bias vs sd

long2 <- temp%>%mutate(bias_sd_ratio=ifelse(sd_surve!=0,mean_surve_bias/sd_surve,0))
plot_ly(long2,x=~t1,y=~t2,z=~abs(bias_sd_ratio)) # bias/sd

mean(abs(long2$bias_sd_ratio)>2) #propotion of ratio magnitude above 2
hist(abs(long2$bias_sd_ratio))

##################
##################
##################
##################
##Investigation on particular points

bias_point_fun<-function(result_list,truth,position_index){
  bias_point<-NULL
  for(j in 1:length(result_list)){

    bias_temp<-(result_list[[j]][[1]][[1]]$survival_curve/result_list[[j]][[1]][[1]]$survival_curve[1]-truth)[position_index]
    N1_l1_temp<-result_list[[j]][[1]][[2]]
    N2_l1_temp<-result_list[[j]][[1]][[3]]
    bias_point<-rbind(bias_point,c(N1_l1_temp,N2_l1_temp,bias_temp,update=F))
  }

  for(j in 1:length(result_list)){

    bias_temp<-(result_list[[j]][[2]][[1]]$survival_curve/result_list[[j]][[2]][[1]]$survival_curve[1]-truth)[position_index]
    N1_l1_temp<-result_list[[j]][[2]][[2]]
    N2_l1_temp<-result_list[[j]][[2]][[3]]
    bias_point<-rbind(bias_point,c(N1_l1_temp,N2_l1_temp,bias_temp,update=T))
  }
  colnames(bias_point)<-c("N1_l1","N2_l1","bias","update")
  bias_point<-as.data.frame(bias_point)
  return(bias_point)
}




#Any point s1 s2
#Estimated d\lambda_A
s1<-0.6
s2<-0.6

#
position_index<-which(round(t1,2)==s1 & round(t2,2)==s2)
bias_point_data_l1<-bias_point_fun(result.l1.1,truth = truth,position_index)
ggplot(bias_point_data_l1,aes(x=as.factor(update),y=bias))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  ggtitle("residual ("%+% s1 %+% "," %+% s2 %+% ") usual 4 good points")

bias_point_data_l1%>%group_by(update)%>%summarise(sd=sd(bias),bias=mean(bias))


# #
# position_index<-which(round(t1,2)==s1 & round(t2,2)==s2)
# bias_point_data_l1<-bias_point_fun(result.l1.2,truth = truth,position_index)
# ggplot(bias_point_data_l1,aes(x=as.factor(update),y=bias))+
#   geom_boxplot()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
#   ggtitle("residual ("%+% s1 %+% "," %+% s2 %+% ") flex t 4 good points")
#
# bias_point_data_l1%>%group_by(update)%>%summarise(sd=sd(bias),bias=mean(bias))
#
# #
# position_index<-which(round(t1,2)==s1 & round(t2,2)==s2)
# bias_point_data_l1<-bias_point_fun(result.l2.1,truth = truth,position_index)
# ggplot(bias_point_data_l1,aes(x=as.factor(update),y=bias))+
#   geom_boxplot()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
#   ggtitle("residual ("%+% s1 %+% "," %+% s2 %+% ") usual 4 good points TailShrinkage")
#
# bias_point_data_l1%>%group_by(update)%>%summarise(sd=sd(bias),bias=mean(bias))
#
#
# #
# position_index<-which(round(t1,2)==s1 & round(t2,2)==s2)
# bias_point_data_l1<-bias_point_fun(result.l2.2,truth = truth,position_index)
# ggplot(bias_point_data_l1,aes(x=as.factor(update),y=bias))+
#   geom_boxplot()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
#   ggtitle("residual ("%+% s1 %+% "," %+% s2 %+% ") flex t 4 good points TailShrinkage")
#
# bias_point_data_l1%>%group_by(update)%>%summarise(sd=sd(bias),bias=mean(bias))
#
#
#
# #######
# #
# position_index<-which(round(t1,2)==s1 & round(t2,2)==s2)
# bias_point_data_l1<-bias_point_fun(result.l3.1,truth = truth,position_index)
# ggplot(bias_point_data_l1,aes(x=as.factor(update),y=bias))+
#   geom_boxplot()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
#   ggtitle("residual ("%+% s1 %+% "," %+% s2 %+% ") usual 4 good points masking")
#
# bias_point_data_l1%>%group_by(update)%>%summarise(sd=sd(bias),bias=mean(bias))
#
#
# #
# position_index<-which(round(t1,2)==s1 & round(t2,2)==s2)
# bias_point_data_l1<-bias_point_fun(result.l3.2,truth = truth,position_index)
# ggplot(bias_point_data_l1,aes(x=as.factor(update),y=bias))+
#   geom_boxplot()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
#   ggtitle("residual ("%+% s1 %+% "," %+% s2 %+% ") flex t 4 good points masking")
#
# bias_point_data_l1%>%group_by(update)%>%summarise(sd=sd(bias),bias=mean(bias))
#
#
#
# temp_N1<-result.l1.4[[1]][[3]][[1]]
# temp_N2<-result.l1.4[[1]][[3]][[3]]
#
#
# ##################
# ##################
# ##################
# ##################










# #Investigation into the fit
# result_temp<-result.l1.10
# sum_surve_initial<-result_temp[[1]][[1]][[1]]$survival_curve/result_temp[[1]][[1]][[1]]$survival_curve[1]
# sum_surve_update<-result_temp[[1]][[2]][[1]]$survival_curve/result_temp[[1]][[2]][[1]]$survival_curve[1]
# for(i in 2: length(result_temp)){
#   sum_surve_initial<-sum_surve_initial+result_temp[[i]][[1]][[1]]$survival_curve/result_temp[[i]][[1]][[1]]$survival_curve[1]
#   sum_surve_update<-sum_surve_update+result_temp[[i]][[2]][[1]]$survival_curve/result_temp[[i]][[2]][[1]]$survival_curve[1]
# }
# mean_surve_initial1<-sum_surve_initial/length(result_temp)
# mean_surve_update1<-sum_surve_update/length(result_temp)
#
# library(plotly)
# plot_ly(x=t1,y=t2,z=abs(mean_surve_update1-truth)-abs(mean_surve_initial1-truth))
#
#
# plot_ly(x=t1,y=t2,z=mean_surve_initial1-truth)
# plot_ly(x=t1,y=t2,z=mean_surve_update1-truth)
#
#
#


# #We can actually mimic the shape
# plot_ly(x=t1,y=t2,z=mean_surve_initial-truth)
# plot_ly(x=t1,y=t2,z=mean_surve_initial-empirical_truth,add=T)
# ####Try to combine the two plots above using colors
#
# plot_ly(x=t1,y=t2,z=empirical_truth-truth)
#
# plot_ly(x=t1,y=t2,z=empirical_truth)
# plot_ly(x=t1,y=t2,z=mean_surve_initial)



# #Investigation into D*
# canonical_sum_fun<-function(result_list){
#   canonical_temp<-NULL
#   for(j in 1:length(result_list)){
#     canonical_gradient_initial<-as.numeric(unlist(lapply(result_list[[j]][[1]][[6]], `[[`, 5)))
#     canonical_gradient_updated<-result_list[[j]][[2]][[6]]
#     sd_initial<-sqrt(mean(canonical_gradient_initial^2)/length(canonical_gradient_initial))
#     mean_initial<-mean(canonical_gradient_initial)
#     sd_update<-sqrt(mean(canonical_gradient_updated^2)/length(canonical_gradient_updated))
#     mean_update<-mean(canonical_gradient_updated)
#     canonical_temp<-rbind(canonical_temp,c(mean_initial,sd_initial,mean_update,sd_update))
#   }
#
#   colnames(canonical_temp)<-c("mean_initial","sd_initial","mean_update","sd_update")
#   canonical_temp<-as.data.frame(canonical_temp)
#   return(canonical_temp)
# }
#
# canonical_data<-canonical_sum_fun(result.l1)
#
# #Investigation into the correlation
#
# temp_result<-result.l1.1
# position_index<-which(round(t1,2)==0.2 & round(t2,2)==0.2)
# bias_point_data_l1<-bias_point_fun(temp_result,truth = truth,position_index)
# temp1<-bias_point_data_l1$bias[bias_point_data_l1$update==0]+truth[position_index]
# temp1.1<-bias_point_data_l1$bias[bias_point_data_l1$update==1]+truth[position_index]
#
# position_index<-which(round(t1,2)==0.2 & round(t2,2)==0.8)
# bias_point_data_l1<-bias_point_fun(temp_result,truth = truth,position_index)
# temp2<-bias_point_data_l1$bias[bias_point_data_l1$update==0]+truth[position_index]
# temp2.1<-bias_point_data_l1$bias[bias_point_data_l1$update==1]+truth[position_index]
#
# position_index<-which(round(t1,2)==0.8 & round(t2,2)==0.2)
# bias_point_data_l1<-bias_point_fun(temp_result,truth = truth,position_index)
# temp3<-bias_point_data_l1$bias[bias_point_data_l1$update==0]+truth[position_index]
# temp3.1<-bias_point_data_l1$bias[bias_point_data_l1$update==1]+truth[position_index]
#
# position_index<-which(round(t1,2)==0.8 & round(t2,2)==0.8)
# bias_point_data_l1<-bias_point_fun(temp_result,truth = truth,position_index)
# temp4<-bias_point_data_l1$bias[bias_point_data_l1$update==0]+truth[position_index]
# temp4.1<-bias_point_data_l1$bias[bias_point_data_l1$update==1]+truth[position_index]
#
# temp_data<-cbind(temp1,temp2,temp3,temp4)
# temp_data.1<-cbind(temp1.1,temp2.1,temp3.1,temp4.1)
#
# cor(temp_data)
# cor(temp_data.1)
#
# #Some clever covariates investigation
# result.l1[[1]][[6]]
#
# hist(unlist(lapply(result.l1_true[[3]][[1]][[6]], `[[`, 5)))
#
# temp<-NULL
# for(i in 1:10){
#   temp<-c(temp,mean((unlist(lapply(result.l1[[i]][[1]][[6]], `[[`, 5)))))
#
# }
# hist(temp)
#
#
#
# temp<-NULL
# for(i in 1:10){
#   temp<-c(temp,unlist(lapply(result.l1[[i]][[1]][[6]], `[[`, 3)))
#
# }
# hist(temp)
#
#
#
#
#
# #For the first points the clever covariates for each O_i and then combine for all O_i
# hist(unlist(lapply(result.l1[[1]][[1]][[6]][[1]], `[[`, 3)))
# hist(unlist(lapply(result.l1[[1]][[1]][[6]][[1]], `[[`, 4)))
#
#
# #PnD* at initial for four points
# mean((unlist(lapply(result.l1[[1]][[1]][[6]][[1]], `[[`, 5))))
# mean((unlist(lapply(result.l1[[1]][[1]][[6]][[2]], `[[`, 5))))
# mean((unlist(lapply(result.l1[[1]][[1]][[6]][[3]], `[[`, 5))))
# mean((unlist(lapply(result.l1[[1]][[1]][[6]][[4]], `[[`, 5))))
#
# mean((unlist(lapply(result.l1[[2]][[1]][[6]][[1]], `[[`, 5))))
# mean((unlist(lapply(result.l1[[2]][[1]][[6]][[2]], `[[`, 5))))
# mean((unlist(lapply(result.l1[[2]][[1]][[6]][[3]], `[[`, 5))))
# mean((unlist(lapply(result.l1[[2]][[1]][[6]][[4]], `[[`, 5))))






#
# #Look into relaxed HAL's clever covariates
# temp.1<-result.l1[[1]][[3]]
# temp.2<-result.l1[[2]][[3]]
#
#
#






# #Some visualization
# ggplot(bias_sd_grid_data_l1,aes(x=as.factor(update),y=N1_l1))+
#   geom_boxplot()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
#   ggtitle("N1_l1")
#
# ggplot(bias_sd_grid_data_l1,aes(x=as.factor(update),y=N2_l1))+
#   geom_boxplot()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
#   ggtitle("N2_l1")
#
