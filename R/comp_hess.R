Hess_e02 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- -1* sum(weight*pr*(1-pr))

  h
}

Hess_emax2 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- -1*sum(weight* pr*(1-pr)*(data$dose/(exp(led50)+data$dose))^2)

  h
}

Hess_led502 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- -1 * sum(weight * (pr*(1-pr)*(emax*data$dose*exp(led50)/(exp(led50)+data$dose)^2)^2+
                            (data$y-pr)*emax*data$dose*exp(led50)*(data$dose-exp(led50))/(exp(led50)+data$dose)^3))

  h
}

Hess_e0emax <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- -1* sum(weight* pr*(1-pr)*data$dose/(exp(led50)+data$dose))

  h
}

Hess_e0led50 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- sum(weight* pr*(1-pr)*data$dose*emax/(exp(led50)+data$dose)^2*exp(led50))

  h
}

Hess_emaxled50 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <-  sum(weight* (pr*(1-pr)*data$dose^2*emax/(exp(led50)+data$dose)^3-(data$y-pr)*data$dose/(exp(led50)+data$dose)^2)*exp(led50))

  h
}


#' @title Compute analytical form of Hessian matrix of Binary Emax model
#' @description Compute Hessian matrix of Binary Emax model
#' @param data A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.
#' @param theta Parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @param weight A vector of working weight from the EM iterations.
#' @return A matrix for Heissian.
#' @details Compute Hessian matrix of Binary Emax model which will be used in finding derivative of Jeffery's prior.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Comp_Hess
#' @export
Comp_Hess <- function(data,theta,weight){
  theta <- as.numeric(theta)
  Hess <- matrix(0,3,3)
  Hess[1,1] <- Hess_e02(data,theta,weight)
  Hess[2,2] <- Hess_emax2(data,theta,weight)
  Hess[3,3] <- Hess_led502(data,theta,weight)
  Hess[1,2] <- Hess_e0emax(data,theta,weight)
  Hess[1,3] <- Hess_e0led50(data,theta,weight)
  Hess[2,3] <- Hess_emaxled50(data,theta,weight)
  Hess[3,1] <- Hess[1,3]
  Hess[3,2] <- Hess[2,3]
  Hess[2,1] <- Hess[1,2]

  Hess
}


#' @title Compute analytical form of expected information matrix of Binary Emax model
#' @description Compute expected information matrix of Binary Emax model
#' @param data A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.
#' @param theta Parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @param weight A vector of working weight from the EM iterations.
#' @return A matrix of expected information
#' @details Compute expected information matrix of Binary Emax model which will be used in finding derivative of Jeffery's prior.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Comp_I
#' @export
Comp_I <- function(data,weight, theta){
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

    I_i <- I_i*weight[i]*pr[i]*(1-pr[i])

    I <- I+I_i
  }

  I
}


