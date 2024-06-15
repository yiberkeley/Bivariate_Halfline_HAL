#Once get the conditional intensity estimated, we want to get the canonical gradient and the estimated survival curve.

Rejection_Sample_DensityBased<-function(coef_N1_initial,basis_list_N1_select,
                     coef_N2_initial,basis_list_N2_select,
                     coef_A1_initial,basis_list_A1_select,
                     coef_A2_initial,basis_list_A2_select,
                     censoring_fit,
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
        if(censoring_fit){
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

        #Get the conditional density without normalizing
        ##indirect approach
        #delt1=0,delt2=0
        scaled_density<-function(t1,t2){
          prob_overall<-1
          if(t1<=t2){

            #Before t1
            ##Create a fine grid from 0,t1, including t1
            if(t1!=0){
              prob_overall<-prob_overall*product_integral_func(start=0,end=t1,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
              prob_overall<-prob_overall*product_integral_func(start=0,end=t1,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
            }
            #At t1
            prob_overall<-prob_overall*hazard_func(t1,0,0,0,basis_list_N1_select,coef_N1,type="N1")

            #Between t1 and t2
            ##Create a fine grid from t1 to t2
            if(t1<t2){
              prob_overall<-prob_overall*product_integral_func(start=t1,end=t2,1,1,t1,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
              #At t2
              prob_overall<-prob_overall*hazard_func(t2,1,1,t1,basis_list_N2_select,coef_N2,type="N2")

            }else{
              #################
              prob_overall<-prob_overall*hazard_func(t2,0,0,0,basis_list_N2_select,coef_N2,type="N2")
              #################
            }



          }else{
            #Before t2
            ##Create a fine grid from 0,t2, including t2
            if(t2!=0){
              prob_overall<-prob_overall*product_integral_func(start=0,end=t2,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
              prob_overall<-prob_overall*product_integral_func(start=0,end=t2,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
            }
            #At t2
            prob_overall<-prob_overall*hazard_func(t2,0,0,0,basis_list_N2_select,coef_N2,type="N2")

            #Between t2 and t1
            ##Create a fine grid from t2 to t1
            prob_overall<-prob_overall*product_integral_func(start=t2,end=t1,1,1,t2,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")

            #At t1
            prob_overall<-prob_overall*hazard_func(t1,1,1,t2,basis_list_N1_select,coef_N1,type="N1")

          }
          return(prob_overall)
        }


        #Use rejection sampling technique
        grid_scale_density<-0.01
        grid_temp<- seq(t_start,1,grid_scale_density)[-1]

        if(half_line_type=="N2"){
            maxDens = max(unlist(lapply(grid_temp,function(x)scaled_density(t_cut,x))))
            size_temp=0
            final_sample<-NULL
            final_sample_weight<-NULL
            while(size_temp<sample_size){
              sampled <- data.frame(proposal = runif(sample_size,t_start,1))
              sampled$targetDensity <-unlist(lapply(sampled$proposal,function(x)scaled_density(t_cut,x)))
              sampled$accepted = ifelse(runif(sample_size,0,1) < sampled$targetDensity / maxDens,
                                        TRUE, FALSE)

              final_sample<-c(final_sample,sampled$proposal[sampled$accepted])
              final_sample_weight<-c(final_sample_weight,sampled$targetDensity[sampled$accepted])
              size_temp<-length(final_sample)
            }

            final_sample<-final_sample[1:sample_size]
            final_sample_weight<-final_sample_weight[1:sample_size]
        }else{
            maxDens = max(unlist(lapply(grid_temp,function(x)scaled_density(x,t_cut))))
            size_temp=0
            final_sample<-NULL
            final_sample_weight<-NULL
            while(size_temp<sample_size){
              sampled <- data.frame(proposal = runif(sample_size,t_start,1))
              sampled$targetDensity <-unlist(lapply(sampled$proposal,function(x)scaled_density(x,t_cut)))
              sampled$accepted = ifelse(runif(sample_size,0,1) < sampled$targetDensity / maxDens,
                                        TRUE, FALSE)

              final_sample<-c(final_sample,sampled$proposal[sampled$accepted])
              final_sample_weight<-c(final_sample_weight,sampled$targetDensity[sampled$accepted])
              size_temp<-length(final_sample)
            }

            final_sample<-final_sample[1:sample_size]
            final_sample_weight<-final_sample_weight[1:sample_size]
        }


        result_inner<-list()
        result_inner[[1]]<-final_sample
        result_inner[[2]]<-final_sample_weight
      return(result_inner)
}


