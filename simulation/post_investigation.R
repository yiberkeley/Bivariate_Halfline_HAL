#In no censoring case, the EIC is 1(t1>s1,t2>s2)-p(s1,s2), So solving the empirical \hat{D}^* is to solve
#the difference between empricial survival curve and updated survival curve.

result.full.indep<-result.full.indep.detail

t1_plot_detail<-t1.detail[seq(1,9000,by=1000)]
t2_plot_detail<-t2.detail[seq(1,9000,by=1000)]

temp_table<-as.data.frame(cbind(t1.detail,t2.detail,index=1:length(t2.detail)))
temp_table<-temp_table%>%filter(t1.detail %in% t1_plot_detail)%>%filter(t2.detail %in% t2_plot_detail)
index<-temp_table$index

seeds<-c(101,500,2000,4321,8688,10113)

for(i in 1:length(seeds)){

    #Get fitted survival
    surve_fit<-result.full.indep[[i]][[2]][[1]]$survival_curve

    print(mean(result.full.indep[[i]][[2]][[6]][[1]]))
    print(mean(result.full.indep[[i]][[2]][[6]][[2]]))
    print(mean(result.full.indep[[i]][[2]][[6]][[3]]))
    print(mean(result.full.indep[[i]][[2]][[6]][[4]]))

    seed<-seeds[i]
    set.seed(seed)
    datgen<-function(n){return(datgen_independent_no_censoring(n))}
    data_full<-datgen(100)
    data_obs<-data_full%>%select(tt1,delt1,tt2,delt2)

    #Get empricial
    t1_plot<- t1_plot_detail
    t2_plot <- t2_plot_detail
    plot_data<-expand.grid(t1=t1_plot,t2=t2_plot)
    empirical<-future_apply(plot_data,1,function(x) mean(data_obs$tt1>x[1] & data_obs$tt2>x[2]))

    # #Get truth
    # set.seed(100)
    # nCores<-detectCores()-2
    # plan(multisession,workers = nCores)
    # data_full_sim<-datgen(1000)
    # truth<-future_apply(plot_data,1,function(x) mean(data_full_sim$t1>x[1] & data_full_sim$t2>x[2]))


    print(summary(empirical-surve_fit[index]))
}

library(plotly)
surve_fit<-result.full.indep[[2]][[2]][[1]]$survival_curve
plot_table<-as.data.frame(cbind(t1=t1.detail,t2=t2.detail,surve=surve_fit))
plot_ly(plot_table,x=~t1,y=~t2,z=~surve)
plot_ly(plot_table%>%filter(t1<=0.8,t2<=0.8),x=~t1,y=~t2,z=~surve)

plot_table2<-plot_table%>%mutate(surve_modified=ifelse((surve-max(plot_table$surve)+1)>=0,(surve-max(plot_table$surve)+1),0))
plot_table2<-plot_table2%>%mutate(truth=(1-t1)*(1-t2))
plot_ly(plot_table2,x=~t1,y=~t2,z=~surve_modified)
plot_ly(plot_table2%>%filter(t1<=0.8,t2<=0.8),x=~t1,y=~t2,z=~surve_modified)
plot_ly(plot_table2,x=~t1,y=~t2,z=~surve_modified-truth)
plot_ly(plot_table2%>%filter(t1<=0.8,t2<=0.8),x=~t1,y=~t2,z=~surve_modified-truth)


seed<-seeds[2]
set.seed(seed)
datgen<-function(n){return(datgen_independent_no_censoring(n))}
data_full<-datgen(100)
data_obs<-data_full%>%select(tt1,delt1,tt2,delt2)

#Get empricial
plot_data3<-as.data.frame(cbind(t1=t1.detail,t2=t2.detail))
empirical<-future_apply(plot_data3,1,function(x) mean(data_obs$tt1>x[1] & data_obs$tt2>x[2]))



plot_ly(x=plot_table2$t1,y=plot_table2$t2,z=plot_table2$surve_modified-empirical)
plot_ly(x=plot_table2$t1,y=plot_table2$t2,z=plot_table2$truth-empirical)



