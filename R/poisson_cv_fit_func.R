
#Poisson cv fit helper function
poisson_cv_fit_func<-function(poisson_data_N3,basis_list_N3,nfolds,penalty.factor_N3,undersmooth=F){
  #CV part
  up<-(-4)
  low<-(-8)
  dist_lambda<-(-0.1)
  lambda_seq<-exp(seq(up,low,dist_lambda))
  stop_cv<-FALSE
  zoom_in_times<-1
  while(stop_cv==FALSE & zoom_in_times <6){
    print("check point: prepare to get cv poisson fit")
    fitness<-cv_fit(poisson_data_N3,lambda_seq,nfolds,penalty.factor_N3,cv_choice="lambda")
    print("check point: cv poisson fit finished")
    num_coef_cv<-fitness$cv_num_coef[which(fitness$cv_loss==min(fitness$cv_loss))]
    if ( sum((unique(num_coef_cv)==sort(fitness$cv_num_coef)[1]) |
         (unique(num_coef_cv)==sort(fitness$cv_num_coef)[2]))>=1 ){
      up_temp<-log(lambda_seq[length(num_coef_cv)])
      up<-up_temp+1/zoom_in_times
      low<-up_temp-1/zoom_in_times
      lambda_seq<-exp(seq(up,low,dist_lambda/zoom_in_times))
      zoom_in_times<-zoom_in_times+1
      print("cv zoom in up")
    }else if (sum((unique(num_coef_cv)==sort(fitness$cv_num_coef,decreasing = T)[1]) |
              (unique(num_coef_cv)==sort(fitness$cv_num_coef,decreasing = T)[2]))>=1 ){
      low_temp<-log(lambda_seq[length(lambda_seq)-length(num_coef_cv)+1])
      up<-low_temp+1
      low<-low_temp-1
      lambda_seq<-exp(seq(up,low,dist_lambda/zoom_in_times))
      zoom_in_times<-zoom_in_times+1
      print("cv zoom in low")
    }else{
      stop_cv<-TRUE
      print("cv finished")
    }
  }

  lambda_chosen_seq=min(exp(fitness$log_lambda[which(fitness$cv_loss==min(fitness$cv_loss))])) #only cv
  cv_selected_cv_risk_N3<-min(fitness$cv_loss)
  cv_selected_cv_number_coef_N3<-unique(fitness$cv_num_coef[which(fitness$cv_loss==min(fitness$cv_loss))])

  fit_N3 <- glmnet(x = as.matrix(poisson_data_N3[,-c(1:9)]),
                   y = poisson_data_N3$jump,
                   family = poisson(),
                   offset = log(poisson_data_N3$timeIntervalLength),
                   thresh=1E-6,maxit=10^6,
                   lambda = lambda_chosen_seq,
                   penalty.factor = penalty.factor_N3,
                   solver="coordinate_descent",
                   weights= poisson_data_N3$weights)

  mu_fit<-predict(fit_N3,
                  newx=as.matrix(poisson_data_N3[,-c(1:9)]),
                  newoffset =log(poisson_data_N3$timeIntervalLength),
                  type = "response")

  cv_selected_fit_risk_N3<-dev_function_1(y=poisson_data_N3$jump,mu=mu_fit,weights=poisson_data_N3$weights)/2

  coef_N3_initial_list<-list()
  basis_list_N3_select_list<-list()
  penalty_N3_select_list<-list()
  for(i in 1:length(lambda_chosen_seq)){
    coef_N3_initial_list[[i]]<-coef(fit_N3)[,i][which(coef(fit_N3)[,i]!=0)]
    basis_list_N3_select_list[[i]]<-basis_list_N3[(which(coef(fit_N3)[,i]!=0)-1)[-1]]
    penalty_N3_select_list[[i]]<-penalty.factor_N3[(which(coef(fit_N3)[,i]!=0)-1)[-1]]
  }

  #Undersmoothing part
  if(undersmooth==T){
    up<-log(lambda_chosen_seq)-0.1
    low<-log(lambda_chosen_seq)-2
    dist_lambda<-(-0.2)
    lambda_seq<-exp(seq(up,low,dist_lambda))
    stop_undersmooth<-FALSE
    zoom_in_times<-1
    while(stop_undersmooth==FALSE & zoom_in_times <4){
      print("check point: prepare to get undersmooth poisson fit")
      fitness<-cv_fit(poisson_data_N3,lambda_seq,nfolds,penalty.factor_N3,cv_choice="lambda")
      print("check point: undersmooth poisson fit finished")
      if ( sum(fitness$cv_num_coef>cv_selected_cv_number_coef_N3*1.5)==0){  #| sum(fitness$cv_loss>2*cv_selected_cv_risk_N3)==0
        up<-log(lambda_seq[length(lambda_seq)])
        low<-up-2
        lambda_seq<-exp(seq(up,low,dist_lambda/zoom_in_times))
        zoom_in_times<-zoom_in_times+1
        print("undersmooth zoom in low")
      }else{
        stop_undersmooth<-TRUE
        print("undersmooth finished")
      }
    }
    lambda_chosen_undersmooth=max(exp(fitness$log_lambda[which( fitness$cv_num_coef>cv_selected_cv_number_coef_N3*1.5)])) #fitness$cv_loss>2*cv_selected_cv_risk_N3 |

    fit_N3_undersmooth <- glmnet(x = as.matrix(poisson_data_N3[,-c(1:9)]),
                     y = poisson_data_N3$jump,
                     family = poisson(),
                     offset = log(poisson_data_N3$timeIntervalLength),
                     thresh=1E-6,maxit=10^6,
                     lambda = lambda_chosen_undersmooth,
                     penalty.factor = penalty.factor_N3,
                     solver="coordinate_descent",
                     weights= poisson_data_N3$weights)

    mu_fit_undersmooth<-predict(fit_N3_undersmooth,
                    newx=as.matrix(poisson_data_N3[,-c(1:9)]),
                    newoffset =log(poisson_data_N3$timeIntervalLength),
                    type = "response")


    coef_N3_initial_list_undersmooth<-list()
    basis_list_N3_select_list_undersmooth<-list()
    penalty_N3_select_list_undersmooth<-list()
    for(i in 1:length(lambda_chosen_undersmooth)){
      coef_N3_initial_list_undersmooth[[i]]<-coef(fit_N3_undersmooth)[,i][which(coef(fit_N3_undersmooth)[,i]!=0)]
      basis_list_N3_select_list_undersmooth[[i]]<-basis_list_N3[(which(coef(fit_N3_undersmooth)[,i]!=0)-1)[-1]]
      penalty_N3_select_list_undersmooth[[i]]<-penalty.factor_N3[(which(coef(fit_N3_undersmooth)[,i]!=0)-1)[-1]]
    }
  }


  result_list<-list()
  result_list[[1]]<-cv_selected_cv_risk_N3
  result_list[[2]]<-cv_selected_cv_number_coef_N3
  result_list[[3]]<-cv_selected_fit_risk_N3
  result_list[[4]]<-coef_N3_initial_list
  result_list[[5]]<-basis_list_N3_select_list
  result_list[[6]]<-penalty_N3_select_list
  if(undersmooth==T){
    result_list[[7]]<-coef_N3_initial_list_undersmooth
    result_list[[8]]<-basis_list_N3_select_list_undersmooth
    result_list[[9]]<-penalty_N3_select_list_undersmooth
  }
  return(result_list)
}
