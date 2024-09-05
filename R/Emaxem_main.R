#' @title Fitting IL method with Emax model and binary response missing data.
#' @description This provides the estimates using IL method.
#' @param data Dataset contain missing response indicated as'NA', including response variable as 'y' and dose variable as 'dose'.
#' @param theta_0 Initial value of Emax model parameters for optimization, Default: NULL
#' @param alpha_0 Initial value of logistic missingness model parameters for optimization., Default: NULL
#' @param mis_form an object of class "formula": a symbolic description of the model to be fitted, Default: as.formula(mis ~ y + dose + x1 + x2)
#' @return list of fitted values:
#' \item{theta}{the final fitted parameters of Emax model}
#' \item{alpha}{the final fitted parameters of logistic missing model}
#' \item{weight}{the final fitted weight for each observation in EM}
#' \item{Q}{the value of Q function for maximizationfor each iteration of EM}
#' \item{K}{the total number of iterations of EM to converge}
#' \item{vcov_theta}{the estimated variance covariance matrix of theta}
#' \item{vcov_alpha}{the estimated variance covariance matrix of alpha}
#' @examples
#' \dontrun{
#' if(interactive()){
#'  theta_true=matrix(c(qlogis(0.1),qlogis(0.8)-qlogis(0.1),log(7.5)),1,3)
#'  colnames(theta_true)<- c('e_0','emax','led_50')
#'  theta_true <- as.data.frame(theta_true)
#'  dose_set <- c(0,7.5,22.5,75,225)
#'  n=355
#'  alpha_true = c(0.5,1,-0.5,0,0) #mis 15 typical
#'  data <-sim_data(theta_true,n,dose_set,alpha_true)
#'  res <- fitEmaxEM(data=data$data,mis_form=as.formula(mis~y+dose) )
#'  }
#' }
#' @seealso
#'  \code{\link[clinDR]{fitEmax}}
#' @rdname fitEmaxEM
#' @export
#' @importFrom clinDR fitEmax
#' @import formula.tools


fitEmaxEM <- function (data=NULL,theta_0=NULL,alpha_0=NULL, mis_form=as.formula(mis~y+dose+x1+x2)){
  options(warn=-1)
  #agument missing data
  data <- Augment_Missing(data)


  #calculate inital Q, theta_0, alpha_0, and weight
  Q <- 0

  weight = rep(1,length(data$y))
  weight[which(data$mis==1&data$y==1)]=sum(data$y==1&data$mis==0)/sum(data$mis==0)
  weight[which(data$mis==1&data$y==0)]=sum(data$y==0&data$mis==0)/sum(data$mis==0)

  fit_mis <- glm(formula=mis_form, data=data,family = binomial)
  alpha_0 <- coef(fit_mis)

  mis_theta<-clinDR::fitEmax(y=data[which(data$mis==0),]$y,dose=data[which(data$mis==0),]$dose,modType=3,binary = TRUE,diagnostics=FALSE)
  if(is.null(mis_theta)){theta_0 <-  startEmax(y=data$y,dose=data$dose,binary = TRUE)}
  else{theta_0 = mis_theta$fit$estimate}
  # theta_0 <-startEmax(data$y,data$dose,binary = TRUE,count = weight)

  theta_0 <- as.data.frame(matrix(c(theta_0[3],theta_0[2],(theta_0[1])),1,3))
  colnames(theta_0)<-c('e_0','emax','led_50')
  theta_m <- theta_0

  theta_res <- comp_theta(data=data,weight = weight,theta=theta_0)
  theta = as.data.frame(t(theta_res$par))
  colnames(theta)<-c('e_0','emax','led_50')
  theta_m <- rbind(theta_m,theta)

  alpha = alpha_0

  Q[2] <- comp_Q(data,theta=theta,alpha=alpha,weight = weight,mis_form=mis_form)

  k <- 2

  while(abs(Q[k]-Q[k-1])>=1e-4){
    # E step
    weight <- comp_weight(data,theta=theta,alpha=alpha,mis_form = mis_form)

    # M step
    form.c <- as.character(mis_form)
    fit_mis <-eval(parse(text = paste("glm(", form.c, ",  weights = weight, data=data,family = binomial)")))
    #estimate alpha via weighted logistic regression
    #fit_mis <- glm(formula = mis_form,  weights = weight, data=data,family = binomial)

    alpha <- coef(fit_mis)

    #estimate alpha via optimization
    # alpha_res <- comp_alpha(data=data,weight = weight,alpha=alpha)
    # if(alpha_res$convergence!=0){
    #   message('alpha does not converge.')
    #   break
    # }
    # alpha <- alpha_res$par

    #estimate theta via optimization

    theta_res <- comp_theta(data=data,weight = weight,theta=theta)
    # if(theta_res$convergence!=0){
    #   message('theta does not converge.')
    #   break
    #   }
    theta <- as.data.frame(t(theta_res$par))
    hessian <- theta_res$hessian
    colnames(theta)<-c('e_0','emax','led_50')
    k <- k+1
    Q[k] <- comp_Q(data,theta=theta,alpha=alpha,weight = weight,mis_form = mis_form)
    theta_m <- rbind(theta_m,theta)
  }



  res <- list(theta=theta,alpha=alpha,weight=weight,Q=Q,K=k,theta_m=theta_m,hessian=hessian,data=data,mis_fit=fit_mis)
  vcov <- comp_vcov(res)
  res$vcov_theta <- vcov$vcov_theta
  res$vcov_alpha <- vcov$vcov_alpha
  return(res)
}




#' @title Fitting bias reduced IL method with Emax model and binary response missing data.
#' @description This provides the estimates using FIL method.
#' @param data Dataset contain missing response indicated as'NA', including response variable as 'y' and dose variable as 'dose'.
#' @param theta_0 Initial value of Emax model parameters for optimization, Default: NULL
#' @param alpha_0 Initial value of logistic missingness model parameters for optimization., Default: NULL
#' @param mis_form an object of class "formula": a symbolic description of the model to be fitted, Default: as.formula(mis ~ y + dose + x1 + x2)
#' @return list of fitted values:
#' \item{theta}{the final fitted parameters of Emax model}
#' \item{alpha}{the final fitted parameters of logistic missing model}
#' \item{weight}{the final fitted weight for each observation in EM}
#' \item{Q}{the value of Q function for maximizationfor each iteration of EM}
#' \item{K}{the total number of iterations of EM to converge}
#' \item{vcov_theta}{the estimated variance covariance matrix of theta}
#' \item{vcov_alpha}{the estimated variance covariance matrix of alpha}
#' @examples
#' \dontrun{
#' if(interactive()){
#'  theta_true=matrix(c(qlogis(0.1),qlogis(0.8)-qlogis(0.1),log(7.5)),1,3)
#'  colnames(theta_true)<- c('e_0','emax','led_50')
#'  theta_true <- as.data.frame(theta_true)
#'  dose_set <- c(0,7.5,22.5,75,225)
#'  n=355
#'  alpha_true = c(0.5,1,-0.5,0,0) #mis 15 typical
#'  data <-sim_data(theta_true,n,dose_set,alpha_true)
#'  res <- fitEmaxEM_firth(data=data$data,mis_form=as.formula(mis~y+dose) )
#'  }
#' }
#' @seealso
#'  \code{\link[clinDR]{fitEmax}}
#' @rdname fitEmaxEM_firth
#' @export
#' @importFrom clinDR fitEmax
#' @import brglm brglm
#' @import formula.tools
fitEmaxEM_firth <- function (data=NULL,theta_0=NULL,alpha_0=NULL,mis_form=as.formula(mis~y+dose+x1+x2)){
  options(warn=-1)
  #agument missing data
  data <- Augment_Missing(data)


  #calculate inital Q, theta_0, alpha_0, and weight
  Q <- 0

  weight = rep(1,length(data$y))

  fit_mis <- brglm::brglm(mis_form, data=data,family = binomial,method = "brglm.fit")
  alpha_0 <- coef(fit_mis)

  mis_theta<-clinDR::fitEmax(y=data[which(data$mis==0),]$y,dose=data[which(data$mis==0),]$dose,modType=3,binary = TRUE,diagnostics=FALSE)
  if(is.null(mis_theta)){theta_0 <-  startEmax(y=data$y,dose=data$dose,binary = TRUE)}
  else{theta_0 = mis_theta$fit$estimate}
  # theta_0 <-startEmax(data$y,data$dose,binary = TRUE,count = weight)

  theta_0 <- as.data.frame(matrix(c(theta_0[3],theta_0[2],(theta_0[1])),1,3))
  colnames(theta_0)<-c('e_0','emax','led_50')
  theta_m <- theta_0
  weight <- comp_weight(data=data, theta = theta_0, alpha = alpha_0,mis_form=mis_form)
  theta_res <- comp_theta_firth(data=data,weight = weight,theta=theta_0)
  theta = as.data.frame(t(theta_res$par))
  colnames(theta)<-c('e_0','emax','led_50')
  theta_m <- rbind(theta_m,theta)

  form.c <- as.character(mis_form)
  fit_mis <-eval(parse(text = paste("brglm::brglm(", form.c, ",  weights = weight, data=data,family = binomial)")))
  #fit_mis = brglm(mis~y+dose+x1+x2, data=data,family = binomial, weights = weight,method = "brglm.fit")
  alpha = coef(fit_mis)

  Q[2] <- comp_Q_firth(data,theta=theta,alpha=alpha,weight = weight,fit_mis=fit_mis)

  k <- 2

  #while(abs(sum((theta_m[k,])^2)-sum((theta_m[k-1,])^2))>=1e-4){
  while(abs(Q[k]-Q[k-1])>=1e-3){
    # E step
    weight <- comp_weight(data,theta=theta,alpha=alpha,mis_form=mis_form)

    # M step

    #estimate alpha via weighted logistic regression
    #fit_mis <- glm(mis~y+dose+x1+x2, data=data,family = binomial,weights=weight)

    form.c <- as.character(mis_form)
    fit_mis <-eval(parse(text = paste("brglm::brglm(", form.c, ",  weights = weight, data=data,family = binomial)")))
    #fit_mis <- brglm(mis~y+dose+x1+x2, data=data,family = binomial, weights = weight,method = "brglm.fit")

    #fit_mis <- brglm(mis~y+dose+x1+x2, data=data,family = binomial,method = "brglm.fit")
    alpha <- coef(fit_mis)

    #estimate alpha via optimization
    # alpha_res <- comp_alpha(data=data,weight = weight,alpha=alpha)
    # if(alpha_res$convergence!=0){
    #   message('alpha does not converge.')
    #   break
    # }
    # alpha <- alpha_res$par

    #estimate theta via optimization
    theta_res <- comp_theta_firth(data=data,weight = weight,theta=theta_0)
    # if(theta_res$convergence!=0){
    #   message('theta does not converge.')
    #   break
    #   }
    theta <- as.data.frame(t(theta_res$par))
    hessian <- theta_res$hessian
    colnames(theta)<-c('e_0','emax','led_50')
    k <- k+1
    Q[k] <- comp_Q_firth(data,theta=theta,alpha=alpha,weight = weight,fit_mis=fit_mis)
    theta_m <- rbind(theta_m,theta)
  }
  weight_j=weight


  ##compute vcov matrix
  res <- list(theta=theta,alpha=alpha,weight=weight,Q=Q,K=k,theta_m=theta_m,hessian=hessian,data=data,mis_fit=fit_mis)
  vcov <- comp_vcov_firth(res)
  res$vcov_theta <- vcov$vcov_theta
  res$vcov_alpha <- vcov$vcov_alpha

  return(res)
}
