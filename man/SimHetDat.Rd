\name{SimHetDat}
\alias{SimHetDat}
\title{ Simulate Heterogeneous Data for Gaussian Graphical Models}
\description{Simulate Heterogeneous data with a band structure, which can be used in \code{GGMM(data,...)} for estimating the structure of the Gaussian graphical network.}

\usage{
SimHetDat(n = 100, p = 200, M = 3, mu = 0.3, type = "band")
}

\arguments{
  \item{ n }{ Number of observations for each group, default of 100.}
  \item{ p }{ Number of covariates for each observation, default of 200.}
  \item{ M }{ Number of latent groups for the simulated dataset choose 2 or 3, default of 3.}
    \item{ mu }{ The mean difference among groups. If \eqn{M=3}, the mean of three groups are \eqn{-mu,0,mu}, respectively. If \eqn{M=2}, the mean of two groups are \eqn{0,mu}, respectively.}
   \item{ type }{\code{type=="band"} which denotes the band structure, with precision matrix
  \deqn{
        C_{i,j}=\left\{\begin{array}{ll}
                     0.5,&\textrm{if $\left| j-i \right|=1, i=2,...,(p-1),$}\\
                     0.25,&\textrm{if $\left| j-i \right|=2, i=3,...,(p-2),$}\\
                     1,&\textrm{if $i=j, i=1,...,p,$}\\
                    0,&\textrm{otherwise.}
                \end{array}\right. }
                }
}
\value{
   \item{ data }{ \emph{n}x\emph{p} Heterogeneous Gaussian distributed data.}
  \item{ A }{  \emph{p}x\emph{p} adjacency matrix used for generating data.}
  \item{label}{The group indices for each observation.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
library(GGMM)
SimHetDat(n = 100, p = 200, M = 3, mu = 0.5, type = "band")
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

}

\keyword{SimHetDat}
