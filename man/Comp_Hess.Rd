% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comp_hess.R
\name{Comp_Hess}
\alias{Comp_Hess}
\title{Compute analytical form of Hessian matrix of Binary Emax model}
\usage{
Comp_Hess(data, theta, weight)
}
\arguments{
\item{data}{A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.}

\item{theta}{Parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).}

\item{weight}{A vector of working weight from the EM iterations.}
}
\value{
A matrix for Heissian.
}
\description{
Compute Hessian matrix of Binary Emax model
}
\details{
Compute Hessian matrix of Binary Emax model which will be used in finding derivative of Jeffery's prior.
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
