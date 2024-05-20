####HAS BEEN CHANGED TO RIDGE TEMPORARILY: NO

# #L1 penalization
# cvrx_poisson_lambda<-function(lambda_seq,x_basis_poisson,offset,weights,y){
#
#   lambda_vals<-lambda_seq
#   p<-dim(x_basis_poisson)[2]
#   n<-dim(x_basis_poisson)[1]
#   # beta_vals <- matrix(0, nrow = p, ncol = length(lambda_vals))
#   # beta_0_vals <- matrix(0, nrow = 1, ncol = length(lambda_vals))
#   X<-x_basis_poisson
#
#   beta_0<-Variable(1)
#   beta <- Variable(p)
#   loss<-sum(exp(offset+X%*% beta+beta_0)*weights*(1-y))-sum((X%*% beta+beta_0)*y*weights)
#
#
#   results <- lapply(lambda_vals,
#                     function (lambda) {
#                       reg <-lambda *p_norm(beta, 1)
#                       obj <- loss + reg
#                       prob <- Problem(Minimize(obj))
#                       result <- solve(prob,solver="ECOS")
#                       result_return<-list()
#                       result_return[[1]]<-result$getValue(beta)
#                       result_return[[2]]<-result$getValue(beta_0)
#                       result_return[[3]]<-result$value
#                       return(result_return)
#                     })
#   return(results)
# }

#L1 constraints with initial offset
cvrx_poisson_l1_offset<-function(l1_norm_seq,x_basis_poisson,offset,weights,y,penalty_factor,time_vec,uniform_ratio,length_clever,coef_initial){
  print("cvrx_poisson_l1_offset_l1_norm_seq" %+% l1_norm_seq)
  p<-dim(x_basis_poisson)[2]
  p<-as.numeric(p)
  n<-dim(x_basis_poisson)[1]
  length_clever<-as.numeric(length_clever)
  # beta_vals <- matrix(0, nrow = p, ncol = length(lambda_vals))
  # beta_0_vals <- matrix(0, nrow = 1, ncol = length(lambda_vals))
  print(p)
  print(length_clever)
  print(p-length_clever)

  X<-x_basis_poisson
  penalty<-penalty_factor
  if((p-length_clever)>=1){
    beta_0<-as.numeric(coef_initial)[1]
    beta_fixed<-as.numeric(coef_initial)[-1]
    epsilon<- Variable(length_clever)

    fixed_part<-X[,1:(p-length_clever)]%*% beta_fixed
    X_var<-X[,(p-length_clever+1):p]
  }else{
    beta_0<-as.numeric(coef_initial)[1]
    epsilon<- Variable(length_clever)
    fixed_part<-0
    X_var<-X[,1:p]
  }

  if(!uniform_ratio){
    loss<-sum(exp(offset+fixed_part+X_var%*%epsilon+beta_0)*weights*(1-y))-sum((fixed_part+X_var%*%epsilon+beta_0)*y*weights)
  }
  else{
    loss<-sum(exp(offset+fixed_part+X_var%*%epsilon+beta_0+log(1/(1-time_vec)))*weights*(1-y))-sum((fixed_part+X_var%*%epsilon+beta_0+log(1/(1-time_vec)))*y*weights)
  }

  results <- lapply(l1_norm_seq,
                    function (l1_norm) {
                      constraint1<-sum(abs(epsilon))<= l1_norm #constraint1<-sum((epsilon)^2)<= l1_norm
                      constraints <-list(constraint1)
                      objective<-Minimize(loss)
                      print("to solve")
                      problem <- Problem(objective, constraints)
                      #browser()
                      result <- CVXR::solve(problem,solver="SCS") #This is crucial ECOS
                      print("solved")
                      result_return<-list()
                      if((p-length_clever)>=1){
                        result_return[[1]]<-c(beta_fixed,result$getValue(epsilon))
                      }else{
                        result_return[[1]]<-c(result$getValue(epsilon))
                      }
                      result_return[[2]]<-beta_0
                      result_return[[3]]<-result$value
                      print("stored")
                      return(result_return)
                    })


  return(results)
}

#L1 constraints
cvrx_poisson_l1<-function(l1_norm_seq,x_basis_poisson,offset,weights,y,penalty_factor,time_vec,uniform_ratio){
  p<-dim(x_basis_poisson)[2]
  n<-dim(x_basis_poisson)[1]
  # beta_vals <- matrix(0, nrow = p, ncol = length(lambda_vals))
  # beta_0_vals <- matrix(0, nrow = 1, ncol = length(lambda_vals))
  X<-x_basis_poisson
  penalty<-penalty_factor

  beta_0<-Variable(1)
  beta <- Variable(p)
  if(!uniform_ratio){
    loss<-sum(exp(offset+X%*% beta+beta_0)*weights*(1-y))-sum((X%*% beta+beta_0)*y*weights)
  }
  else{
    loss<-sum(exp(offset+X%*% beta+beta_0+log(1/(1-time_vec)))*weights*(1-y))-sum((X%*% beta+beta_0+log(1/(1-time_vec)))*y*weights)
  }

  results <- lapply(l1_norm_seq,
                      function (l1_norm) {
                      constraint1<-sum(abs(beta)*penalty)<= l1_norm#constraint1<-sum((beta)^2*penalty)<= l1_norm
                      constraints <-list(constraint1)
                      objective<-Minimize(loss)
                      print("to solve")
                      problem <- Problem(objective, constraints)
                      #browser()
                      result <- CVXR::solve(problem,solver="SCS") #This is crucial ECOS
                      print("solved")
                      result_return<-list()
                      result_return[[1]]<-result$getValue(beta)
                      result_return[[2]]<-result$getValue(beta_0)
                      result_return[[3]]<-result$value
                      print("stored")
                      return(result_return)
                    })


  return(results)
}



cv_fit_cvxr<-function(poisson_data,l1_norm_seq,nfolds,penalty_factor,uniform_ratio=F,position_basis=NULL,additional_offset=NULL,
                      offset_initial=NULL,length_clever=NULL,coef_initial=NULL){

  nfolds<-nfolds
  repeated_data_all<-poisson_data

  #num_cores <- detectCores()-2
  #registerDoParallel(cores = num_cores)
  if(is.null(position_basis)){
    penalty_factor<-penalty_factor
  }else{
    penalty_factor<-penalty_factor[position_basis]
  }

  if(is.null(additional_offset)){
    additional_offset<-rep(0,dim(repeated_data_all)[[1]])
  }else{
    additional_offset<-additional_offset
  }


  print("check point: cv fit prepared")

  cv_results<-NULL
  for(i in 1:nfolds){
    train<-repeated_data_all[repeated_data_all$fold!=i,]
    test<-repeated_data_all[repeated_data_all$fold==i,]

    additional_offset_train<-additional_offset[repeated_data_all$fold!=i]
    additional_offset_test<-additional_offset[repeated_data_all$fold==i]


    if(is.null(position_basis)){
      x_basis_poisson <- as.matrix(train[,-c(1:9)])
      x_basis_poisson_test <- as.matrix(test[,-c(1:9)])
    }else
    {
      x_basis_poisson <- as.matrix(train[,-c(1:9)])[,position_basis]
      x_basis_poisson_test <- as.matrix(test[,-c(1:9)])[,position_basis]
    }

    #Training
    offset<-as.vector(log(train$timeIntervalLength)+additional_offset_train)
    weights<- as.vector(train$weights)
    y<-as.vector(train$jump)
    time_vec<-as.vector(train$t)

    if(is.null(offset_initial)){
      fit<-cvrx_poisson_l1(l1_norm_seq,x_basis_poisson,offset,weights,y,penalty_factor,time_vec,uniform_ratio)
    }else{
      fit<-cvrx_poisson_l1_offset(l1_norm_seq,x_basis_poisson,offset,weights,y,penalty_factor,time_vec,uniform_ratio,length_clever,coef_initial)
    }


    #Validation
    loss_results<-NULL
    for(j in 1:length(l1_norm_seq)) {


      offset_test<-log(test$timeIntervalLength)+additional_offset_test
      weights_test<- test$weights
      y_test<-test$jump
      X_test<-x_basis_poisson_test
      time_vec_test<-as.vector(test$t)

      beta<-fit[[j]][[1]]
      beta_0<-fit[[j]][[2]]
      if(!uniform_ratio){
        loss_valid<-sum(exp(offset_test+X_test%*% beta+beta_0)*weights_test*(1-y_test))-sum((X_test%*% beta+beta_0)*y_test*weights_test)
      }else{
        loss_valid<-sum(exp(offset_test+X_test%*% beta+beta_0+log(1/(1-time_vec_test)))*weights_test*(1-y_test))-sum((X_test%*% beta+beta_0+log(1/(1-time_vec_test)))*y_test*weights_test)
      }


      cv_num_coef_ij <- sum(abs(beta)>10^-4)
      ##########
      print("dev1 within cv" %+% loss_valid)
      #PnL
      cv_loss_ij <-loss_valid
      cv_l1_norm_ij<-sum(abs(beta)[abs(beta)>10^-4])
      loss_results<-rbind(loss_results,c(l1_norm_seq[j],cv_num_coef_ij, cv_loss_ij,cv_l1_norm_ij,i))
    }

    cv_results<-rbind(cv_results,loss_results)
    #return(loss_results)
  }

  colnames(cv_results)<-c('l1_norm','cv_num_coef','cv_loss','l1','round')
  cv_results<-as.data.frame(cv_results)

  #library(dplyr)

    cv_results_final<-cv_results %>% group_by(l1_norm)%>%summarise(cv_num_coef=mean(cv_num_coef),
                                                                      cv_loss=mean(cv_loss),
                                                                      cv_l1_norm=mean(l1),
                                                                      l1_norm=mean(l1_norm))
  return(cv_results_final)
}

