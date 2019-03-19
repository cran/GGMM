\name{BRGM}
\alias{BRGM}
\title{ Learning gene regulatory networks for breast cancer.}
\description{Gaussian Graphical Mixture Models for learning gene regulatory network with multiple subtypes of breat cancer dataset.}

\usage{
BRGM(data, M=3, alpha1 = 0.05, alpha2 = 0.02, alpha3 = 0.2, iteration = 30, warm = 20)
}

\arguments{
  \item{ data }{ A \emph{\eqn{n}}x\emph{\eqn{p}} matrix of breast cancer expression data.}
    \item{ M }{ The number of heterogeneous groups, default of 3 based on the BIC scores.}
    \item{ alpha1  }{ The significance level of correlation screening in the \eqn{\psi}-learning algorithm, see R package \pkg{equSA} for detail. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set, which reduces the risk of missing important variables in 
 the conditioning set. In general, including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient, the default value is 0.05.}
  \item{ alpha2  }{ The significance level of \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix, see \pkg{equSA}, the default value is 0.02.}
\item{ alpha3  }{ The significance level of integrative \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix of GGMM method, the default value is 0.2.}
\item{ iteration  }{ The number of total iterations, the default value is 30.}
\item{ warm  }{ The number of burn-in iterations, the default value is 20.}
}

\value{
  \item{ Adj }{  \emph{\eqn{p}}x\emph{\eqn{p}} Estimated adjacency matrix for network construction.}
  \item{ label }{ The estimated group indices for each observation.}
  \item{ BIC }{ The BIC scores for determining the number of groups \eqn{M}.}
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
 \dontshow{
suppressWarnings(RNGversion("3.5.0"))
set.seed(100)
library(GGMM)
data(breast)
Gra <- BRGM(breast[1:60,1:5], M=2,alpha1 = 0.8, alpha2 = 0.9, alpha3 = 0.9, iteration = 1, warm = 0)
}
\donttest{

library(GGMM)
library(huge)
data(breast)
## long time ##
M <- 3
Gra <- BRGM(breast, M, iteration = 30, warm = 20)
## plot gene regulatory network by our estimated adjacency matrix.
huge.plot(Gra$Adj)

}



         
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2018). An Imputation Regularized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{BRGM}
