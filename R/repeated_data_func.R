
#A function to create repeated data

repeated_data_func<-function(data_obs,foldid,perc_sample,monitoring_time,censor_data_include=T,masking_option,weight_option){
        ##Repeated Data with modeling information for N1,A1
        repeated_data_N1<-NULL
        repeated_data_A1<-NULL
        for(i in 1:dim(data_obs)[1]){
          #if tt1 happens first
          if(data_obs$tt1[i]<=data_obs$tt2[i]){

            #Monitoring time up till tt1
            time_point<-monitoring_time[which(monitoring_time<=data_obs$tt1[i])]

            #If there are a lot of time points, then do a random sampling of the monitoring time
            if(length(time_point)>20){
              time_point<-time_point[c(1,sample(2:(length(time_point)-1), floor(perc_sample*(length(time_point)-2)),
                                                replace=FALSE),length(time_point))]
              ############################
              #add the point 0.7 if exceed
              if(data_obs$tt1[i]>0.7){
                time_point<-c(time_point,0.7)
              }
              time_point<-sort(unique(time_point))
            }

            #Get the length of time where conditional intensity will stay the same
            time_interval<-c(diff(time_point),exp(-10)) #Here, the perturbation is needed for the poisson fit from 0 to exp(-10)

            #Count is named after poisson regression, but it basically captures whether the failure happened in the time interval or not
            #If it is a failure outcome, we let the count to be 1 at last time point.
            count<-c(rep(0,length(time_point)-1),1)
            if(data_obs$delt1[i]==1){
              repeated_data_N1_indiv<-data.frame(obs=i,fold=foldid[i],monitoringTime=time_point,timeIntervalLength=time_interval,
                                                 jump=count,tt2=data_obs$tt2[i],delt2=data_obs$delt2[i])
              repeated_data_N1<-rbind(repeated_data_N1, repeated_data_N1_indiv)
              #It may happen that tt1=0
              if(length(time_interval)>1){
                repeated_data_A1_indiv<-repeated_data_N1_indiv[-dim(repeated_data_N1_indiv)[1],]
                repeated_data_A1<-rbind(repeated_data_A1, repeated_data_A1_indiv)
              }
            }
            else{
              repeated_data_A1_indiv<-data.frame(obs=i,fold=foldid[i],monitoringTime=time_point,timeIntervalLength=time_interval,
                                                 jump=count,tt2=data_obs$tt2[i],delt2=data_obs$delt2[i])
              repeated_data_A1<-rbind(repeated_data_A1, repeated_data_A1_indiv)
              #It may happen that tt1=0
              if(length(time_interval)>1){
                repeated_data_N1_indiv<-repeated_data_A1_indiv[-dim(repeated_data_A1_indiv)[1],]
                repeated_data_N1<-rbind(repeated_data_N1, repeated_data_N1_indiv)
              }
            }
          }
          else{

            #If tt1 happens later than tt2
            #if further tt1 is a failure time
            time_point_2<-monitoring_time[which(monitoring_time<=data_obs$tt1[i])]

            #Subsampling the monitoring time before tt2 and after tt2 seperately
            change_point_index<-which(time_point_2==data_obs$tt2[i])

            if(change_point_index>20){
              temp_sample<-c(1,sample(2:(change_point_index-1), floor(perc_sample*(change_point_index-2)), replace=FALSE),change_point_index)
            }else{
              temp_sample<-1:change_point_index
            }

            if(length(time_point_2)-change_point_index>20){
              temp_sample_2<-c(sample((change_point_index+1): (length(time_point_2)-1),
                                      floor(perc_sample*( length(time_point_2)-1-(change_point_index+1))),
                                      replace=FALSE),length(time_point_2))
            }else{
              temp_sample_2<-(change_point_index+1):length(time_point_2)
            }

            time_point_2<-time_point_2[c(temp_sample,temp_sample_2)]
            ############################
            #add the point 0.7 if exceed
            if(data_obs$tt1[i]>0.7){
              time_point_2<-c(time_point_2,0.7)
            }
            time_point_2<-sort(unique(time_point_2))


            time_interval_2<-c(diff(time_point_2),exp(-10))
            count_2<-c(rep(0,length(time_point_2)-1),1)
            if(data_obs$delt1[i]==1){
              repeated_data_N1_indiv<-data.frame(obs=i,fold=foldid[i],monitoringTime=time_point_2,timeIntervalLength=time_interval_2,
                                                 jump=count_2,tt2=data_obs$tt2[i],delt2=data_obs$delt2[i])
              repeated_data_N1<-rbind(repeated_data_N1, repeated_data_N1_indiv)

              repeated_data_A1_indiv<-repeated_data_N1_indiv[-dim(repeated_data_N1_indiv)[1],]
              repeated_data_A1<-rbind(repeated_data_A1, repeated_data_A1_indiv)
            }
            else{
              repeated_data_A1_indiv<-data.frame(obs=i,fold=foldid[i],monitoringTime=time_point_2,timeIntervalLength=time_interval_2,
                                                 jump=count_2,tt2=data_obs$tt2[i],delt2=data_obs$delt2[i])
              repeated_data_A1<-rbind(repeated_data_A1, repeated_data_A1_indiv)

              repeated_data_N1_indiv<-repeated_data_A1_indiv[-dim(repeated_data_A1_indiv)[1],]
              repeated_data_N1<-rbind(repeated_data_N1, repeated_data_N1_indiv)
            }
          }
        }


        ##Repeated Data with modeling information for N2,A2
        repeated_data_N2<-NULL
        repeated_data_A2<-NULL
        for(i in 1:dim(data_obs)[1]){
          #if tt2 happens first
          if(data_obs$tt2[i]<=data_obs$tt1[i]){

            #Monitoring time up till tt2
            time_point<-monitoring_time[which(monitoring_time<=data_obs$tt2[i])]

            #If there are a lot of time points, then do a random sampling of the monitoring time
            if(length(time_point)>20){
              time_point<-time_point[c(1,sample(2:(length(time_point)-1), floor(perc_sample*(length(time_point)-2)),
                                                replace=FALSE),length(time_point))]
              ############################
              #add the point 0.7 if exceed
              if(data_obs$tt2[i]>0.7){
                time_point<-c(time_point,0.7)
              }
              time_point<-sort(unique(time_point))
            }

            #Get the length of time where conditional intensity will stay the same
            time_interval<-c(diff(time_point),exp(-10)) #Here, the perturbation is needed for the poisson fit from 0 to exp(-10)

            #Count is named after poisson regression, but it basically captures whether the failure happened in the time interval or not
            #If it is a failure outcome, we let the count to be 1 at last time point.
            count<-c(rep(0,length(time_point)-1),1)
            if(data_obs$delt2[i]==1){
              repeated_data_N2_indiv<-data.frame(obs=i,fold=foldid[i],monitoringTime=time_point,timeIntervalLength=time_interval,
                                                 jump=count,tt1=data_obs$tt1[i],delt1=data_obs$delt1[i])
              repeated_data_N2<-rbind(repeated_data_N2, repeated_data_N2_indiv)
              #It may happen that tt2=0
              if(length(time_interval)>1){
                repeated_data_A2_indiv<-repeated_data_N2_indiv[-dim(repeated_data_N2_indiv)[1],]
                repeated_data_A2<-rbind(repeated_data_A2, repeated_data_A2_indiv)
              }
            }
            else{
              repeated_data_A2_indiv<-data.frame(obs=i,fold=foldid[i],monitoringTime=time_point,timeIntervalLength=time_interval,
                                                 jump=count,tt1=data_obs$tt1[i],delt1=data_obs$delt1[i])
              repeated_data_A2<-rbind(repeated_data_A2, repeated_data_A2_indiv)
              #It may happen that tt2=0
              if(length(time_interval)>1){
                repeated_data_N2_indiv<-repeated_data_A2_indiv[-dim(repeated_data_A2_indiv)[1],]
                repeated_data_N2<-rbind(repeated_data_N2, repeated_data_N2_indiv)
              }
            }
          }
          else{

            #If tt2 happens later than tt1
            time_point_2<-monitoring_time[which(monitoring_time<=data_obs$tt2[i])]

            #Subsampling the monitoring time before tt1 and after tt1 seperately
            change_point_index<-which(time_point_2==data_obs$tt1[i])

            if(change_point_index>20){
              temp_sample<-c(1,sample(2:(change_point_index-1), floor(perc_sample*(change_point_index-2)), replace=FALSE),change_point_index)
            }else{
              temp_sample<-1:change_point_index
            }

            if(length(time_point_2)-change_point_index>20){
              temp_sample_2<-c(sample((change_point_index+1): (length(time_point_2)-1),
                                      floor(perc_sample*( length(time_point_2)-1-(change_point_index+1))),
                                      replace=FALSE),length(time_point_2))
            }else{
              temp_sample_2<-(change_point_index+1):length(time_point_2)
            }

            time_point_2<-time_point_2[c(temp_sample,temp_sample_2)]
            ############################
            #add the point 0.7 if exceed
            if(data_obs$tt2[i]>0.7){
              time_point_2<-c(time_point_2,0.7)
            }
            time_point_2<-sort(unique(time_point_2))
            time_interval_2<-c(diff(time_point_2),exp(-10))
            count_2<-c(rep(0,length(time_point_2)-1),1)
            if(data_obs$delt2[i]==1){
              repeated_data_N2_indiv<-data.frame(obs=i,fold=foldid[i],monitoringTime=time_point_2,timeIntervalLength=time_interval_2,
                                                 jump=count_2,tt1=data_obs$tt1[i],delt1=data_obs$delt1[i])
              repeated_data_N2<-rbind(repeated_data_N2, repeated_data_N2_indiv)

              repeated_data_A2_indiv<-repeated_data_N2_indiv[-dim(repeated_data_N2_indiv)[1],]
              repeated_data_A2<-rbind(repeated_data_A2, repeated_data_A2_indiv)
            }
            else{
              repeated_data_A2_indiv<-data.frame(obs=i,fold=foldid[i],monitoringTime=time_point_2,timeIntervalLength=time_interval_2,
                                                 jump=count_2,tt1=data_obs$tt1[i],delt1=data_obs$delt1[i])
              repeated_data_A2<-rbind(repeated_data_A2, repeated_data_A2_indiv)

              repeated_data_N2_indiv<-repeated_data_A2_indiv[-dim(repeated_data_A2_indiv)[1],]
              repeated_data_N2<-rbind(repeated_data_N2, repeated_data_N2_indiv)
            }
          }
        }





        ##Repeated Data Step 2
        #Repeated Data Step 2 for N1,A1
        #I(tt2<t)
        augmentaion_part_1_N1<-as.numeric(repeated_data_N1$tt2<repeated_data_N1$monitoringTime)
        augmentaion_part_1_A1<-as.numeric(repeated_data_A1$tt2<repeated_data_A1$monitoringTime)
        #I(tt2<t)*delt2
        augmentaion_part_2_N1<-augmentaion_part_1_N1*repeated_data_N1$delt2
        augmentaion_part_2_A1<-augmentaion_part_1_A1*repeated_data_A1$delt2
        #I(tt2<t)tt2
        augmentaion_part_3_N1<-augmentaion_part_1_N1*repeated_data_N1$tt2
        augmentaion_part_3_A1<-augmentaion_part_1_A1*repeated_data_A1$tt2

        repeated_data_N1_2<-cbind(repeated_data_N1,
                                  t=repeated_data_N1$monitoringTime,
                                  tt2_less_t=augmentaion_part_1_N1,
                                  tt2_less_t_delt2=augmentaion_part_2_N1,
                                  tt2_less_t_tt2=augmentaion_part_3_N1)
        repeated_data_N1_2<-repeated_data_N1_2[,-c(3,6,7)]

        repeated_data_A1_2<-cbind(repeated_data_A1,
                                  t=repeated_data_A1$monitoringTime,
                                  tt2_less_t=augmentaion_part_1_A1,
                                  tt2_less_t_delt2=augmentaion_part_2_A1,
                                  tt2_less_t_tt2=augmentaion_part_3_A1)
        repeated_data_A1_2<-repeated_data_A1_2[,-c(3,6,7)]


        #Repeated Data Step 2 for N2,A2
        #I(tt1<t)
        augmentaion_part_1_N2<-as.numeric(repeated_data_N2$tt1<repeated_data_N2$monitoringTime)
        augmentaion_part_1_A2<-as.numeric(repeated_data_A2$tt1<repeated_data_A2$monitoringTime)
        #I(tt1<t)*delt1
        augmentaion_part_2_N2<-augmentaion_part_1_N2*repeated_data_N2$delt1
        augmentaion_part_2_A2<-augmentaion_part_1_A2*repeated_data_A2$delt1
        #I(tt1<t)tt1
        augmentaion_part_3_N2<-augmentaion_part_1_N2*repeated_data_N2$tt1
        augmentaion_part_3_A2<-augmentaion_part_1_A2*repeated_data_A2$tt1

        repeated_data_N2_2<-cbind(repeated_data_N2,
                                  t=repeated_data_N2$monitoringTime,
                                  tt1_less_t=augmentaion_part_1_N2,
                                  tt1_less_t_delt1=augmentaion_part_2_N2,
                                  tt1_less_t_tt1=augmentaion_part_3_N2)
        repeated_data_N2_2<-repeated_data_N2_2[,-c(3,6,7)]

        repeated_data_A2_2<-cbind(repeated_data_A2,
                                  t=repeated_data_A2$monitoringTime,
                                  tt1_less_t=augmentaion_part_1_A2,
                                  tt1_less_t_delt1=augmentaion_part_2_A2,
                                  tt1_less_t_tt1=augmentaion_part_3_A2)
        repeated_data_A2_2<-repeated_data_A2_2[,-c(3,6,7)]



        if (masking_option=="masking"){
          repeated_data_N1_2<-repeated_data_N1_2%>%filter(t<=0.7)
          repeated_data_N2_2<-repeated_data_N2_2%>%filter(t<=0.7)
          repeated_data_A1_2<-repeated_data_A1_2%>%filter(t<=0.7)
          repeated_data_A2_2<-repeated_data_A2_2%>%filter(t<=0.7)
        }

        if (weight_option=="tail_shrinkage"){
          repeated_data_N1_2<-repeated_data_N1_2%>%mutate(weights=ifelse(t<=0.7,1,0.1))
          repeated_data_N2_2<-repeated_data_N2_2%>%mutate(weights=ifelse(t<=0.7,1,0.1))
          repeated_data_A1_2<-repeated_data_A1_2%>%mutate(weights=ifelse(t<=0.7,1,0.1))
          repeated_data_A2_2<-repeated_data_A2_2%>%mutate(weights=ifelse(t<=0.7,1,0.1))
        }else if(weight_option=="no_tail_shrinkage"){
          repeated_data_N1_2<-repeated_data_N1_2%>%mutate(weights=1)
          repeated_data_N2_2<-repeated_data_N2_2%>%mutate(weights=1)
          repeated_data_A1_2<-repeated_data_A1_2%>%mutate(weights=1)
          repeated_data_A2_2<-repeated_data_A2_2%>%mutate(weights=1)
          }

        result_list<-list()

        if(censor_data_include){
          result_list[[1]]<-repeated_data_N1_2
          result_list[[2]]<-repeated_data_N2_2
          result_list[[3]]<-repeated_data_A1_2
          result_list[[4]]<-repeated_data_A2_2
        }else{
          result_list[[1]]<-repeated_data_N1_2
          result_list[[2]]<-repeated_data_N2_2
        }

        return(result_list)

}

