# Expectation-Maximization Based Emax Model Package

In this document, we illustrate the main features of the `ememax` R package through examples. Additional information on the statistical methodology and computational details are provided in the accompanying documentation and research articles.

## Cite the package

The package applies methods introduced in the [paper](https://arxiv.org/abs/2410.00259):

Zhang, J., Pradhan, V. and Zhao, Y., 2024. Robust Emax Model Fitting: Addressing Nonignorable Missing Binary Outcome in Dose-Response Analysis. arXiv preprint arXiv:2410.00259.

## Install

Open the R console and run the following command to install the package from source:

```r
install.packages("devtools") # When you have not installed devtools package
devtools::install_github("Celaeno1017/ememax")
```

## Tutorial

First, load the R package.

```r
library(ememax)
```

To illustrate the main features of the R package `ememax`, let's first generate some data. We have built in a few functions directly into the R package for this purpose.

```r
theta_true=matrix(c(qlogis(0.1),qlogis(0.8)-qlogis(0.1),log(7.5)),1,3) #true parmaeter of emax model
 colnames(theta_true)<- c('e_0','emax','led_50')
 theta_true <- as.data.frame(theta_true)
 dose_set <- c(0,7.5,22.5,75,225) #doseage definition
 n=355 #total number of sample size. The sample will be evenly allocated.
 alpha_true = c(0.5,1,-0.5,0,0) #mis 15 typical
 data <-sim_data(theta_true,n,dose_set,alpha_true)
```

To fit the emEmax model with Firth type correction, we use the fitEmaxEM_firth function which implements the proposed methodology.

```r
res <- fitEmaxEM_firth(data=data$data,mis_form=as.formula(mis~y+dose) )
```
Key parameters include:
- `mis_form`: The pre-defined logistic model for missingness. If y is included as a covariate, that means one considers the missingness is non-ignorable.

The result will contain the following values:
- `theta`:	the final fitted parameters of Emax model

- `alpha`:	the final fitted parameters of the logistic missing model

- `weight`: the final fitted weight for each observation in EM

- `Q`:	the value of Q function for maximizationfor each iteration of EM

- `K`: the total number of iterations of EM to converge

- `vcov_theta`: the estimated variance covariance matrix of theta

- `vcov_alpha`:	the estimated variance covariance matrix of alpha
