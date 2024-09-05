#' @title Calculate the variance covariance matrix of estimated parameters by EmaxEM
#' @description This provides the estimated VCOV matrix of parameters using IL method by EmaxEM.
#' @param em.emax.fit an object for result from EmaxEM
#' @return A list of two variance covariance matrices, one for Emax model parameter, one for logistic missing model parameter.
#' @details Internal function for variance covariance estimation.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[MASS]{ginv}}
#'  \code{\link[boot]{inv.logit}}
#' @rdname comp_vcov
#' @export
#' @importFrom MASS ginv
#' @importFrom boot inv.logit

comp_vcov <- function(em.emax.fit){
  theta = as.numeric(em.emax.fit$theta)
  weight = em.emax.fit$weight
  hessian = em.emax.fit$hessian
  data = em.emax.fit$data

  hessian = Comp_Hess(data=data,theta=theta,weight = weight)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]
  ######vcov of theta

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  #pr <- exp(emax*data$dose/(ed50+data$dose)+e0)/(1+exp(emax*data$dose/(ed50+data$dose)+e0))
  s_e0 <- (data$y-pr)
  s_emax <- (data$y-pr)*data$dose/(exp(led50)+data$dose)
  s_led50 <- (data$y-pr)*(-1)*data$dose*emax*exp(led50)/(exp(led50)+data$dose)^2

  S <- cbind(s_e0,s_emax,s_led50)

  # second_term <- matrix(0,3,3)
  # for (i in 1:length(data$y)){
  #   SS <- S[i,]%*%t(S[i,])
  #   second_term <- second_term+SS*weight[i]
  # }
  second_term <- t(S)%*%diag(weight)%*%S

  third_term <- matrix(0,3,3)
  for (i in unique(data$id)){
    Q <- c(sum(weight[which(data$id==i)]*s_e0[which(data$id==i)]),
           sum(weight[which(data$id==i)]*s_emax[which(data$id==i)]),
           sum(weight[which(data$id==i)]*s_led50[which(data$id==i)]))
    third_term <- third_term+ Q%*%t(Q)
  }

  #Q <- c(sum(weight*s_e0),sum(weight*s_emax),sum(weight*s_led50))
  Information <-  -hessian - second_term + third_term

  vcov <- MASS::ginv(Information)

  ####vcov of alpha
  mis_fit <- em.emax.fit$mis_fit
  Infm <- solve(vcov(mis_fit))
  current.parameter = coef(mis_fit)
  full.X <- model.matrix(mis_fit$formula,data=data)
  #full.X <- cbind(rep(1,length(data$y)),as.matrix(cbind(data$y,data$dose,data$x1,data$x2)))
  mu <- boot::inv.logit(full.X %*% current.parameter)
  w <-mis_fit$weights
  full.y <- data$y
  full.r <- data$mis
  q <- NULL
  W <- diag(w)
  Ws <- sqrt(W)

  p1 <- dim(full.X)[2]
  n.full <- dim(full.X)[1]

  for (j in 1:p1) {
    q[j] <- sum(w * (full.y - mu) *
                  full.X[, j])
  }
  s <- matrix(0, nrow = n.full, ncol = p1)
  for (i in 1:n.full) {
    for (j in 1:p1) {
      s[i, j] <- (full.y[i] - mu[i] ) * full.X[i, j]
    }
  }
  second.term <- 0
  for (i in 1:n.full) {
    second.term <- second.term + w[i] * (s[i, ] %*% t(s[i,
    ]))
  }
  Information <- Infm - (second.term - q %*% t(q))
  vcov_mis <- MASS::ginv(Information)

  return(list(vcov_theta=vcov,vcov_alpha=vcov_mis))


}



comp_Q_e0 <- function(data, weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)

  Q <- matrix(0,3,3)
  Q[2,3] <- sum(-1*weight^2* pr*(1-pr)*data$dose/(exp(led50)+data$dose)^2*exp(led50))
  Q[3,3] <- sum(weight^2* pr*(1-pr)*(data$dose-exp(led50))*emax*data$dose*exp(led50)/(exp(led50)+data$dose)^3)
  Q[3,2] <- Q[2,3]

  Q
}

comp_Q_emax <- function(data,weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)

  Q <- matrix(0,3,3)
  Q[2,3] <- sum(-1*weight^2* pr*(1-pr)*data$dose^2/(exp(led50)+data$dose)^3*exp(led50))
  Q[3,3] <- sum(weight^2* pr*(1-pr)*(data$dose-exp(led50))*emax*data$dose^2*exp(led50)/(exp(led50)+data$dose)^4)
  Q[3,2] <- Q[2,3]

  Q
}

comp_Q_led50 <- function(data,weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)

  Q <- matrix(0,3,3)
  Q[2,3] <- sum(weight^2* pr*(1-pr)*data$dose^2*emax/(exp(led50)+data$dose)^4*exp(led50)^2)
  Q[3,3] <- sum(-1*weight^2* pr*(1-pr)*(data$dose-exp(led50))*emax^2*data$dose^2*exp(led50)^2/(exp(led50)+data$dose)^5)
  Q[3,2] <- Q[2,3]

  Q
}

comp_P_e0 <- function(data, weight, theta){

  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)


  I <- matrix(0,3,3)
  for(i in 1:length(data$y)){
    I_i <- matrix(0,3,3)
    I_i[1,] <- c(1,data$dose[i]/(exp(led50)+data$dose[i]),-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)
    I_i[2,] <- I_i[1,]*data$dose[i]/(exp(led50)+data$dose[i])
    I_i[3,] <- I_i[2,]*exp(led50)*(-emax)/(exp(led50)+data$dose[i])

    I_i <- I_i*weight[i]^2*(pr[i]*(1-pr[i])^3-pr[i]^3*(1-pr[i]))

    I <- I+I_i
  }

  I

}

comp_P_emax <- function(data, weight, theta){

  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)


  I <- matrix(0,3,3)
  for(i in 1:length(data$y)){
    I_i <- matrix(0,3,3)
    I_i[1,] <- c(1,data$dose[i]/(exp(led50)+data$dose[i]),-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)
    I_i[2,] <- I_i[1,]*data$dose[i]/(exp(led50)+data$dose[i])
    I_i[3,] <- I_i[2,]*exp(led50)*(-emax)/(exp(led50)+data$dose[i])

    I_i <- I_i*weight[i]^2*(pr[i]*(1-pr[i])^3-pr[i]^3*(1-pr[i]))*data$dose[i]/(exp(led50)+data$dose[i])

    I <- I+I_i
  }

  I

}

comp_P_led50 <- function(data, weight, theta){

  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)


  I <- matrix(0,3,3)
  for(i in 1:length(data$y)){
    I_i <- matrix(0,3,3)
    I_i[1,] <- c(1,data$dose[i]/(exp(led50)+data$dose[i]),-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)
    I_i[2,] <- I_i[1,]*data$dose[i]/(exp(led50)+data$dose[i])
    I_i[3,] <- I_i[2,]*exp(led50)*(-emax)/(exp(led50)+data$dose[i])

    I_i <- I_i*weight[i]^2*(pr[i]*(1-pr[i])^3-pr[i]^3*(1-pr[i]))*(-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)

    I <- I+I_i
  }

  I

}

Score_e0 <- function(data, weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  q1 <- sum(weight*(data$y-pr))
  Q_e0 <- comp_Q_e0(data, weight, theta)
  P_e0 <- comp_P_e0(data, weight, theta)
  I <- Comp_I(data, weight, theta)
  I_inv <- MASS::ginv(I)

  U <- q1+sum(diag(I_inv%*%(Q_e0+P_e0)))/2

  U
}

Score_emax <- function(data, weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  q2 <- sum(weight*(data$y-pr)*data$dose/(exp(led50)+data$dose))
  Q_emax <- comp_Q_emax(data, weight, theta)
  P_emax <- comp_P_emax(data, weight, theta)
  I <- Comp_I(data, weight, theta)
  I_inv <- MASS::ginv(I)

  U <- q2+sum(diag(I_inv%*%(Q_emax+P_emax)))/2

  U
}

Score_led50 <- function(data, weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  q3 <- sum(weight*(data$y-pr)*(-1)*data$dose*emax*exp(led50)/(exp(led50)+data$dose)^2)
  Q_led50 <- comp_Q_led50(data, weight, theta)
  P_led50 <- comp_P_led50(data, weight, theta)
  I <- Comp_I(data, weight, theta)
  I_inv <- MASS::ginv(I)

  U <- q3+sum(diag(I_inv%*%(Q_led50+P_led50)))/2

  U
}


#' @title Calculate the variance covariance matrix of estimated parameters by EmaxEM_firth
#' @description This provides the estimated VCOV matrix of parameters using FIL method by EmaxEM_firth.
#' @param em.emax.fit an object for result from EmaxEM_firth
#' @return A list of two variance covariance matrices, one for Emax model parameter, one for logistic missing model parameter.
#' @details Internal function for variance covariance estimation.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[MASS]{ginv}}
#'  \code{\link[boot]{inv.logit}}
#' @rdname comp_vcov_firth
#' @export
#' @importFrom MASS ginv
#' @importFrom boot inv.logit
comp_vcov_firth <- function(em.emax.fit){
  theta = as.numeric(em.emax.fit$theta)
  weight = em.emax.fit$weight
  hessian = em.emax.fit$hessian
  data = em.emax.fit$data


  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  #pr <- exp(emax*data$dose/(ed50+data$dose)+e0)/(1+exp(emax*data$dose/(ed50+data$dose)+e0))

  ##compute s part
  weight_0 = rep(1,length(data$y))
  #weight_0 = weight
  Q_e0 <- comp_Q_e0(data, weight_0, theta)
  Q_emax <- comp_Q_emax(data, weight_0, theta)
  Q_led50 <- comp_Q_led50(data, weight_0, theta)
  P_e0 <- comp_P_e0(data, weight_0, theta)
  P_emax <- comp_P_emax(data, weight_0, theta)
  P_led50 <- comp_P_led50(data, weight_0, theta)
  I <- Comp_I(data, weight_0, theta)
  I_inv <- MASS::ginv(I)

  s_e0 <- (data$y-pr) +sum(diag(I_inv%*%(P_e0+Q_e0)))/(sum(weight)*2)
  s_emax <- (data$y-pr)*data$dose/(exp(led50)+data$dose) +sum(diag(I_inv%*%(P_emax+Q_emax)))/(sum(weight)*2)
  s_led50 <- (data$y-pr)*(-1)*data$dose*emax*exp(led50)/(exp(led50)+data$dose)^2 +sum(diag(I_inv%*%(P_led50+Q_led50)))/(sum(weight)*2)



  S <- cbind(s_e0,s_emax,s_led50)

  second_term <- t(S)%*%diag(weight)%*%S

  ##compute remainder part
  third_term <- matrix(0,3,3)
  for (i in unique(data$id)){
    Q <- c(sum(weight[which(data$id==i)]*s_e0[which(data$id==i)]),
           sum(weight[which(data$id==i)]*s_emax[which(data$id==i)]),
           sum(weight[which(data$id==i)]*s_led50[which(data$id==i)]))
    third_term <- third_term+ Q%*%t(Q)
  }


  # Q <- c(Score_e0(data, weight, theta),Score_emax(data, weight, theta),
  #                 Score_led50(data, weight, theta))
  Information <- -hessian - second_term + third_term

  vcov <- MASS::ginv(Information)
  ####vcov of alpha
  mis_fit <- em.emax.fit$mis_fit
  Infm <- solve(vcov(mis_fit))
  current.parameter = coef(mis_fit)
  full.X <- model.matrix(mis_fit$formula,data=data)
  #full.X <- cbind(rep(1,length(data$y)),as.matrix(cbind(data$y,data$dose,data$x1,data$x2)))
  mu <- boot::inv.logit(full.X %*% current.parameter)
  w <-mis_fit$weights
  full.y <- data$y
  full.r <- data$mis
  q <- NULL
  W <- diag(w)
  Ws <- sqrt(W)
  H <- Ws %*% full.X %*% chol2inv((chol(t(full.X) %*% W %*%
                                          full.X))) %*% t(full.X) %*% Ws
  p1 <- dim(full.X)[2]
  n.full <- dim(full.X)[1]

  for (j in 1:p1) {
    q[j] <- sum(w * (full.y - mu + H[j, j] * (0.5 - mu)) *
                  full.X[, j])
  }
  s <- matrix(0, nrow = n.full, ncol = p1)
  for (i in 1:n.full) {
    for (j in 1:p1) {
      s[i, j] <- (full.y[i] - mu[i] + H[j, j] * (0.5 -
                                                   mu[i])) * full.X[i, j]
    }
  }
  second.term <- 0
  for (i in 1:n.full) {
    second.term <- second.term + w[i] * (s[i, ] %*% t(s[i,
    ]))
  }
  Information <- Infm - (second.term - q %*% t(q))
  vcov_mis <- MASS::ginv(Information)

  return(list(vcov_theta=vcov,vcov_alpha=vcov_mis))


}


