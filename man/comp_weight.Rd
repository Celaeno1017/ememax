% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comp_parameters.R
\name{comp_weight}
\alias{comp_weight}
\title{Estimation of working weight in EM algorithm iteration.}
\usage{
comp_weight(data, theta, alpha, mis_form)
}
\arguments{
\item{data}{A complete dataset without missingness, including response variable as 'y' and dose variable as 'dose'.}

\item{theta}{Parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).}

\item{alpha}{Parameters of logisitic missing indicator model.}

\item{mis_form}{an object of class "formula": a symbolic description of the model to be fitted.}
}
\value{
A vector of weights for each observation in the data.
}
\description{
Internal function for estimating the observation weights for each iteration of EM algorithm.
}
\details{
Calculating the working weights of EM algorithm for each iteration.
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
