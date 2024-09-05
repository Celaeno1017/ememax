#' @title Augment missing response with observed information
#' @description Augment missing response with 0 and 1, and remain all other variables the same.
#' @param data Dataset contain missing response indicated as'NA', including response variable as 'y' and dose variable as 'dose'.
#' @return A complete dataset with augmentation.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Augment_Missing
#' @export

Augment_Missing <- function (data){
  ind.NA <- which(is.na(data$y))
  data$mis <- rep(0,length(data$y))
  data$mis[ind.NA] <- 1


  data <- rbind(data,data[ind.NA,])
  data$y[ind.NA]<- 1
  data$y[which(is.na(data$y))] <- 0

  return(data)
}


#' @title Log likelihood estimation of binary Emax model
#' @description Estimate Log likelihood of given parameters with data for binary Emax model.
#' @param data A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.
#' @param theta Parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @return A vector of log likelihood of each observation.
#' @details Internal function of calculating the maximization function of EM algorithm.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname log_Emax_i
#' @export

log_Emax_i <- function(data,theta){

  log_Emax <- data$y*(theta$e_0+theta$emax*data$dose/(exp(theta$led_50)+data$dose))-
    log(1+exp(theta$e_0+theta$emax*data$dose/(exp(theta$led_50)+data$dose)))
  log_Emax <- plogis((-1)^(data$y+1)*(theta$e_0+theta$emax*data$dose/(exp(theta$led_50)+data$dose)),log.p = TRUE)
  return(log_Emax)
}

#' @title Log likelihood estimation of logisitic missing indicator model
#' @description  Estimate Log likelihood of given parameters with data for logisitic missing indicator model.
#' @param data A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.
#' @param alpha Parameters of logisitic missing indicator model.
#' @param mis_form an object of class "formula": a symbolic description of the model to be fitted.
#' @return A vector of log likelihood of each observation.
#' @details Internal function of calculating the maximization function of EM algorithm.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname log_missing_i
#' @export


log_missing_i <- function(data,alpha,mis_form){
  x <- model.matrix(mis_form,data=data)
  # x <- cbind(rep(1,length(data$y)),as.matrix(cbind(data$y,data$dose,data$x1,data$x2)))
  alpha <- as.matrix(alpha)
  log_missing <- (x%*%alpha)*data$mis-
    log(1+exp(x%*%alpha))
  log_missing <- plogis((-1)^(data$mis+1)*(x%*%alpha),log.p = TRUE)

  return(log_missing)
}


#' @title Maximization function estimation of EM algorithm with defined weight.
#' @description Estimate Maximization function with given parameters, weight, and data for the EM Emax model.
#' @param data A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.
#' @param theta Parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @param alpha Parameters of logisitic missing indicator model.
#' @param weight A vector of working weight from the EM iterations.
#' @param mis_form an object of class "formula": a symbolic description of the model to be fitted.
#' @return A value of function estimation
#' @details Calculating the maximization function of EM algorithm for each iteration.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname comp_Q
#' @export

comp_Q <- function(data,theta,alpha,weight,mis_form){
  log_Emax <- log_Emax_i(data=data,theta=theta)
  log_missing <- log_missing_i(data=data,alpha=alpha,mis_form=mis_form)


  Q = sum(weight*(log_Emax+log_missing))
  return(Q)
}


#' @title Estimation of working weight in EM algorithm iteration.
#' @description Internal function for estimating the observation weights for each iteration of EM algorithm.
#' @param data A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.
#' @param theta Parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @param alpha Parameters of logisitic missing indicator model.
#' @param mis_form an object of class "formula": a symbolic description of the model to be fitted.
#' @return A vector of weights for each observation in the data.
#' @details Calculating the working weights of EM algorithm for each iteration.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname comp_weight
#' @export

comp_weight <- function(data,theta,alpha,mis_form){
  weight <- rep(1,length(data$y))
  log_Emax <- log_Emax_i(data=data,theta=theta)
  log_missing <-log_missing_i(data=data,alpha=alpha,mis_form = mis_form)

  ll <- log_Emax + log_missing

  mis_ind <- which(data$mis==1)
  for (i in mis_ind){
    weight[i] <- exp(ll[i]-log(sum(exp(ll[which(data$id==data$id[i])]))))
  }

  return(weight)

}


#' @title Estimation of emax parameters in EM algorithm iteration.
#' @description Calls Newton-Raphson optimizer MaxNR, for Emax model.
#' @param data A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.
#' @param weight A vector of working weight from the EM iterations.
#' @param theta Initial value of parameters for optimization. The order of the variables is (E0,log(ED50),Emax).
#' @return A vector of parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @details Fits the Emax model with defined working weights using MaxNR.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname comp_theta
#' @export
#' @importFrom maxLik maxNR
comp_theta <- function(data,weight,theta){
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  qr <- function(data,weight,theta){
    e0 <- theta[1]
    emax <- theta[2]
    #led50 <-  ifelse(theta[3]>log(1.2*max(data$dose)),log(1.2*max(data$dose)),theta[3])
    led50 <- theta[3]
    dose <- data$dose
    y <- data$y

    b <- -emax  #emax
    c <- 1/exp(led50) #ed50 related
    lambda <- 1
    a <- e0 - b #e0/emax
    evec <- a + b/(1 + (c * dose)^lambda)

    ll1 <- 0
    ll0 <- 0
    if (length(y == 1) > 0) {
      esub <- evec[y == 1]
      countsub <- weight[y == 1]
      ll1 <- sum(countsub * plogis(esub, log.p = TRUE))
    }
    if (length(y == 0) > 0) {
      esub <- evec[y == 0]
      countsub <- weight[y == 0]
      ll0 <- sum(countsub * plogis(-esub, log.p = TRUE))
    }
    return((ll1 + ll0))
  }

  grr <- function(data,weight,theta){
    e0 <- theta[1]
    emax <- theta[2]
    led50 <- ifelse(theta[3]>log(1.2*max(data$dose)),log(1.2*max(data$dose)),theta[3])

    pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
    #pr <- exp(emax*data$dose/(ed50+data$dose)+e0)/(1+exp(emax*data$dose/(ed50+data$dose)+e0))
    q1 <- sum(weight*(data$y-pr))
    q2 <- sum(weight*(data$y-pr)*data$dose/(exp(led50)+data$dose))
    q3 <- sum(weight*(data$y-pr)*(-1)*data$dose*emax*exp(led50)/(exp(led50)+data$dose)^2)

    c(q1,q2,q3)
  }

  # fit <- nlm(f = qr, p = as.numeric(c(e0,emax,ed50)), hessian = TRUE,
  #            data=data,weight=weight)

  #res<-spg(par=as.numeric(c(e0,emax,ed50)),fn=qr, lower=c(-Inf,-Inf,10^-1), control = list(maxit = 1e5),data=data,weight=weight)
  # res <- constrOptim(theta = as.numeric(c(e0,emax,led50)), f=qr, grad=grr,
  #                    ui=rbind(c(-1,1,0),  # the y-x > 0
  #                             c(0,0,1),   # the z > 0.001*maxdose
  #                             c(0,0,-1)), # the z<1.2*maxdose
  #                    ci=c(0,log(0.001*max(data$dose)),-log(1.2*max(data$dose))),data=data,weight=weight,hessian = TRUE)
  #
  A <- matrix(c(-1,1,0),1,3,byrow = TRUE)
  B <- c(0)

  res <- maxNR(fn=qr, grad=grr,start =as.numeric(c(e0,emax,led50)),data=data,weight=weight )
  return(list(par=res$estimate,hessian=res$hessian))
}



#' @title Estimation of emax parameters in Jeffery's prior penalized IL algorithm iteration.
#' @description Calls Newton-Raphson optimizer MaxNR, for Emax model.
#' @param data A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.
#' @param weight A vector of working weight from the EM iterations.
#' @param theta Initial value of parameters for optimization. The order of the variables is (E0,log(ED50),Emax).
#' @return A vector of parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @details Fits the Emax model with defined working weights using MaxNR with bias reduction.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname comp_theta_firth
#' @export
#' @importFrom maxLik maxNR
comp_theta_firth <- function(data,weight,theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]



  qr <- function(data,weight,theta){
    theta <- as.numeric(theta)
    #inout Hess

    Hess <- Comp_Hess(data,theta,weight)
    #I <- Comp_I(data,weight,theta)
    #logdetI <- ifelse(is.na(log(det(-Hess))),log(det(I)),log(det(-Hess)))
    logdetI <- ifelse(is.na(log(det(-Hess))),0,log(det(-Hess)))

    # I <- Comp_I(data,theta,weight)
    # logdetI <- ifelse(is.na(log(det(I))),0,log(det(I)))
    e0 <- theta[1]
    emax <- theta[2]
    led50 <- theta[3]
    dose <- data$dose
    y <- data$y

    b <- -emax  #emax
    c <- 1/exp(led50) #ed50 related
    lambda <- 1
    a <- e0 - b #e0/emax
    evec <- a + b/(1 + (c * dose)^lambda)

    ll1 <- 0
    ll0 <- 0
    if (length(y == 1) > 0) {
      esub <- evec[y == 1]
      countsub <- weight[y == 1]
      ll1 <- sum(countsub * (plogis(esub, log.p = TRUE)))
    }
    if (length(y == 0) > 0) {
      esub <- evec[y == 0]
      countsub <- weight[y == 0]
      ll0 <- sum(countsub * (plogis(-esub, log.p = TRUE)))
    }
    return((ll1 + ll0+logdetI/2))
  }



  # fit <- nlm(f = qr, p = as.numeric(c(e0,emax,ed50)), hessian = TRUE,
  #            data=data,weight=weight)

  #res<-spg(par=as.numeric(c(e0,emax,ed50)),fn=qr, lower=c(-Inf,-Inf,10^-1), control = list(maxit = 1e5),data=data,weight=weight)
  # res <- constrOptim(theta = as.numeric(c(e0,emax,led50)), f=qr,
  #                    grad= NULL,
  #                    ui=rbind(c(-1,1,0),  # the y-x > 0
  #                             c(0,0,1),   # the z > 0.001*maxdose
  #                             c(0,0,-1)), # the z<1.2*maxdose
  #                    ci=c(0,log(0.001*max(data$dose)),-log(1.2*max(data$dose))),data=data,weight=weight)
  #
  res <- maxNR(fn=qr, grad=NULL,start =as.numeric(c(e0,emax,led50)),data=data,weight=weight )

  #H<-hessian(qr,x= as.numeric(res$par),data=data,weight=weight)
  # return(res)
  return(list(par=res$estimate,hessian=res$hessian))
}


#' @title Maximization function estimation of bias reduced EM algorithm with defined weight.
#' @description Estimate Maximization function with given parameters, weight, and data for the bias reduced EM Emax model.
#' @param data A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.
#' @param theta Parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @param alpha Parameters of logisitic missing indicator model.
#' @param weight A vector of working weight from the EM iterations.
#' @param fit_mis an object of class "formula": a symbolic description of the model to be fitted.
#' @return A value of function estimation
#' @details Calculating the maximization function of bias reduced EM algorithm for each iteration.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname comp_Q_firth
#' @export
comp_Q_firth <- function(data,theta,alpha,weight,fit_mis){

  log_Emax <- log_Emax_i(data=data,theta=theta)
  log_missing <- log_missing_i(data=data,alpha=alpha,mis_form =fit_mis$formula)

  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  #input Hess

  Hess <- Comp_Hess(data,theta,weight)
  #I <- Comp_I(data,weight,theta)
  # logdetI <- ifelse(is.na(log(det(-Hess))),log(det(I)),log(det(-Hess)))
  logdetI <- ifelse(is.na(log(det(-Hess))),0,log(det(-Hess)))
  if(is.na(log(det(-Hess)))){message('det less than 0!')}
  # x <- cbind(rep(1,length(data$y)),as.matrix(cbind(data$y,data$dose,data$x1,data$x2)))
  # alpha <- as.matrix(alpha)
  # pr_mis <- plogis(x%*%alpha)
  # Hess <- -1*t(x)%*%diag(as.numeric(weight*pr_mis*(1-pr_mis)))%*%x
  # logdetI_mis <- log(det(-Hess))

  Q_mis <- fit_mis$deviance/(-2)
  Q = sum(weight*(log_Emax))+logdetI/2+Q_mis
  return(Q)
}



