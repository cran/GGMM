\name{GGMM-package}
\alias{GGMM-package}
\docType{package}
\title{

Gaussian Graphical Mixture Models


}
\description{

The Gaussian graphical model is a widely used tool for learning gene regulatory networks with high-dimensional gene expression data. For many real problems, the data are heterogeneous, which may contain some subgroups or come from different resources. This package provide a Gaussian Graphical Mixture Model (GGMM) for the heterogeneous data. }
\details{
\tabular{rl}{
Package: \tab GGMM\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2019-03-17\cr
License: \tab GPL-2\cr
}



This package illustrates the use of the Gaussian Graphical Mixture Model in two parts:
\bold{}

The first part is to apply the GGMM to estimate network structures using high-dimensional heterogeneous data with a simulated dataset \code{SimHetDat(n,p,...)} and our proposed method \code{GGMM(data,...)}.
\bold{}

The second part is to apply the GGMM to learn a real data example \code{BRGM(breast,...)}, i.e. to learn a common gene regulatory network with heterogeneous gene expression data of breast cancer. The real data example are from The Cancer Genome Atlas (TCGA) with code \code{data(breast)}.

}
\author{
Bochao Jia, Faming Liang
Maintainer: Bochao Jia<jbc409@ufl.edu>
}
\references{


Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.<doi:10.1080/01621459.2015.1012391>

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.<doi:10.1093/biomet/asn036>

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2018). An Imputation Regularized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B. <arXiv:1802.02251>

Jia, B., Xu, S., Xiao, G., Lamba, V., Liang, F. (2017) Inference of Genetic Networks from Next Generation Sequencing Data. Biometrics.

Jia, B. and Liang, F. (2018). Learning Gene Regulatory Networks with High-Dimensional Heterogeneous Data. Accept by ICSA Springer Book. <arXiv:1805.02547>

}
\keyword{ package }

\examples{
 \dontshow{
set.seed(10)
library(GGMM)
result <- SimHetDat(n = 30, p = 10, M = 3, mu = 0, type = "band")
Est <- GGMM(result$data, result$A, alpha1 = 0.5, alpha2 = 0.5, alpha3 = 0.5, M = 3, iteration = 1, warm = 0)
}
\donttest{
library(GGMM)
library(huge)
result <- SimHetDat(n = 100, p = 200, M = 2, mu = 0.5, type = "band")
Graph <- GGMM(result$data, result$A, M = 2, iteration = 30, warm = 20)
      }
}

