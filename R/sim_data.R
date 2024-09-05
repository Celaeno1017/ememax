sim_complete_data <- function(theta,n,dose_set){
  dose <- rep(dose_set,each=round(n/5))
  pi <- plogis(theta$e_0+theta$emax*dose/(exp(theta$led_50)+dose))
  y <- rbinom(n,1,pi)
  x <- as.data.frame(matrix(rnorm(2*n,mean = 1),n,2))
  colnames(x)<-c('x1','x2')
  id <- seq(1,length(y),1)
  data <- cbind(id,y,dose,x)
  return(data)
}

#' @title Simulate dataset for testing Emaxem and Emaxem_firth
#' @description FUNCTION_DESCRIPTION
#' @param theta True parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @param n number of observations.
#' @param dose_set A vector indicate the dose set for the dose-response relationship.
#' @param alpha True parameters of logisitic missing model.
#' @return A list of two datasets. One is with missingness on repsonse, and the other is the full complete data.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  theta_true=matrix(c(qlogis(0.1),qlogis(0.8)-qlogis(0.1),log(7.5)),1,3)
#' colnames(theta_true)<- c('e_0','emax','led_50')
#' theta_true <- as.data.frame(theta_true)
#' dose_set <- c(0,7.5,22.5,75,225)
#' n=355
#' alpha_true = c(0.5,1,-0.5,0,0) #mis 15 typical

#' data <-sim_data(theta_true,n,dose_set,alpha_true)

#'  }
#' }
#' @rdname sim_data
#' @export
sim_data <- function(theta,n,dose_set,alpha){
  comp_data <- sim_complete_data(theta=theta,n=n,dose_set=dose_set)
  comp_data <- as.data.frame(comp_data)
  x <- cbind(rep(1,dim(comp_data)[1]),comp_data$y,comp_data$dose,comp_data$x1,comp_data$x2)
  pi <- plogis(x%*%alpha)
  r <- rbinom(n,1,pi)
  data <- as.data.frame(comp_data)
  data$y[r==1] <- NA
  if (sum((data$y[data$dose==0&r==0]))==0){data$y[sample(which(data$dose==0&r==0),1)]=1}
  return(list(data=data,complete_data=comp_data))
}
