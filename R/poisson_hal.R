
##Fit poisson regression to get the hal fit for the conditional intensities
#Due to the glmnet poisson doesn't allow 0 offset because it uses deviance, we do a perturb and counter strategy
#Modify the dev function to equal to 2*P_nL
#In order to get obj_function as −loglik/nobs+λ∗penalty


dev_function <- function(y, mu, weights, family) {

  if(family$family == "poisson"){
    #2*sum(mu)-2*sum(y*mu) is the 2*[first sum over j<k(i) term of the mean]
    #2*sum((log(mu)+3)*y) is the 2*[sum of f(li,ci) over i]

    #result<-2*sum(mu)-2*sum(y*mu)-2*sum((log(mu)+3)*y) #no weights
    #result<-result/sum(y)

    result<-2*sum(mu*weights)-2*sum(y*mu*weights)-2*sum((log(mu)+10)*y*weights)
    result<-result/sum(y*weights)
  }else
  {result<-sum(family$dev.resids(y, mu, weights))}
  return(result)
}

tmpfun <- get("dev_function",
              envir = asNamespace("glmnet"))
environment(dev_function) <- environment(tmpfun)
assignInNamespace("dev_function",
                  dev_function, ns = "glmnet")



#Built-in cv doesn't work due to no group, create our own cv proceduer

dev_function_1 <- function(y, mu,weights) {

   #result<-2*sum(mu)-2*sum(y*mu)-2*sum((log(mu)+3)*y) #no weights
  #result<-result/sum(y)

  result<-2*sum(mu*weights)-2*sum(y*mu*weights)-2*sum((log(mu)+10)*y*weights)
  result<-result/sum(y*weights)
  return(result)
}

cv_fit<-function(poisson_data,lambda_seq,nfolds,penalty_factor,position_basis=NULL,additional_offset=NULL,cv_choice="lambda"){

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

  # ncores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
  # registerDoParallel(ncores)

  print("check point: cv fit prepared")

  #doFuture::registerDoFuture()

  #cv_results<-foreach(i = 1:nfolds, .combine = rbind) %dopar% {
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

    print("check point: fit for " %+% i)
    fit <- glmnet(x = x_basis_poisson,
                  y = train$jump, family = poisson(),
                  offset = log(train$timeIntervalLength)+additional_offset_train,
                  lambda = lambda_seq,
                  thresh=1E-6,maxit=10^6,
                  penalty.factor = penalty_factor,
                  solver="coordinate_descent",
                  weights = train$weights)

    mu_test<-predict(fit,newx=x_basis_poisson_test,newoffset =log(test$timeIntervalLength)+additional_offset_test,type = "response")

    # loss_results <- foreach(j = 1:dim(mu_test)[2], .combine = rbind) %do% {
    #   cv_num_coef_ij <- sum(!coef(fit, s = lambda_seq[j])[-1] == 0)
    #   loss_valid<-dev_function_1(y=test$jump,mu=mu_test[,j])
    #   #PnL
    #   cv_loss_ij <-loss_valid/2
    #   cv_l1_norm_ij<-sum(abs(coef(fit, s = lambda_seq[j])[-1]))
    #   return(c(log(lambda_seq[j]),cv_num_coef_ij, cv_loss_ij,cv_l1_norm_ij,i))
    # }
    loss_results<-NULL
    for(j in 1:dim(mu_test)[2]) {
      cv_num_coef_ij <- sum(!coef(fit, s = lambda_seq[j])[-1] == 0)
      loss_valid<-dev_function_1(y=test$jump,mu=mu_test[,j],weights=test$weights)
      #PnL
      cv_loss_ij <-loss_valid/2
      cv_l1_norm_ij<-sum(abs(coef(fit, s = lambda_seq[j])[-1]))
      loss_results<-rbind(loss_results,c(log(lambda_seq[j]),cv_num_coef_ij, cv_loss_ij,cv_l1_norm_ij,i))
    }

    cv_results<-rbind(cv_results,loss_results)
    #return(loss_results)
  }

  colnames(cv_results)<-c('log_lambda','cv_num_coef','cv_loss','l1','round')
  cv_results<-as.data.frame(cv_results)

  #library(dplyr)

  if(cv_choice=="lambda"){
    cv_results_final<-cv_results %>% group_by(log_lambda)%>%summarise(cv_num_coef=mean(cv_num_coef),
                                                                      cv_loss=mean(cv_loss),
                                                                      cv_l1_norm=mean(l1),
                                                                      log_lambda1=format(round(mean(log_lambda),10),10))
  }else{
    number_to_gridedNum<-function(grid_start,grid_dist,number){
      return(floor((number-grid_start)/grid_dist)*grid_dist+grid_start)
    }

    grid_start<-min(cv_results$l1)
    grid_dist<-median(diff(sort(cv_results$l1)))*nfolds
    cv_results_temp<-cv_results%>%filter(round==1)
    gridedL1<-sapply(cv_results_temp$l1,function(x){number_to_gridedNum(grid_start,grid_dist,x)})
    cv_results_temp$l1<-gridedL1
    temp<-cv_results_temp%>%group_by(l1)%>%summarise_all(mean)

    for(i in 2:nfolds){
      cv_results_temp<-cv_results%>%filter(round==i)
      gridedL1<-sapply(cv_results_temp$l1,function(x){number_to_gridedNum(grid_start,grid_dist,x)})
      cv_results_temp$l1<-gridedL1
      cv_results_temp_final<-cv_results_temp%>%group_by(l1)%>%summarise_all(mean)

      temp<-inner_join(temp,cv_results_temp_final,by="l1")
    }

    temp_final<-temp %>%mutate(cv_loss= rowMeans(select(., starts_with("cv_loss")), na.rm = T),
                               log_lambda= rowMeans(select(., starts_with("log_lambda")), na.rm = T),
                               cv_num_coef= rowMeans(select(., starts_with("cv_num_coef")), na.rm = T))%>%
                                  select(l1,log_lambda,cv_num_coef,cv_loss)


    cv_results_final<-cbind(temp_final,grid_dist=grid_dist)
  }
  return(cv_results_final)
}

