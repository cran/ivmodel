\name{ivmodel}
\alias{ivmodel}
\title{Fitting Instrumental Variables (IV) Models}
\description{
\code{ivmodel} fits an instrumental variables (IV) model with one endogenous variable and a continuous outcome. It carries out several IV regressions, diagnostics, and tests associated this IV model. It is robust to most data formats, including factor and character data, and can handle very large IV models efficiently.
}
\usage{
ivmodel(Y, D, Z, X, intercept = TRUE, beta0 = 0, alpha = 0.05, 
        k = c(0, 1), heteroSE = FALSE, deltarange = NULL)
}
\arguments{
  \item{Y}{A numeric vector of outcomes}
  \item{D}{A vector of endogenous variables}
  \item{Z}{A matrix or data frame of instruments}
  \item{X}{A matrix or data frame of (exogenous) covariates}
  \item{intercept}{Should the intercept be included? Default is TRUE}
  \item{beta0}{Null value \eqn{\beta_0} for testing null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}. Default is $0$}
  \item{alpha}{The significance level for hypothesis testing. Default is 0.05}
  \item{k}{A numeric vector of k values for k-class estimation. Default is 0 (OLS) and 1 (TSLS)}
  \item{heteroSE}{Should heteroscedastic-robust standard errors be used? Default is FALSE}
  \item{deltarange}{Range of \eqn{\delta} for sensitivity analysis with the Anderson-Rubin test}
}
\details{
Let \eqn{Y_i}, \eqn{D_i}, \eqn{X_{i.}}, and \eqn{Z_{i.}} represent \eqn{i}th individual's outcome, endogenous variable, \emph{p} dimensional exogenous covariates, and \emph{L} dimensional instruments, respectively; note that the intercept can be considered as a vector of ones and a part of the exogenous covariates \eqn{X_{i.}}. \code{ivmodel} assumes the following IV model
\deqn{ Y_i = X_{i.} \alpha + D_i \beta + \epsilon_i, E(\epsilon_i | X_{i.}, Z_{i.}) = 0}
and produces statistics for \eqn{\beta}. In particular, \code{ivmodel} computes the OLS, TSLS, k-class, limited information maximum likelihood (LIML), and Fuller-k  (Fuller 1977) estimates of \eqn{\beta} using \code{kClassESt}, \code{LIMLEst}, and code{FullerEst}. Also, \code{ivmodel} computes confidence intervals and hypothesis tests of the type \eqn{H_0: \beta = \beta_0} versus \eqn{H_0: \beta \neq \beta_0} for the said estimators as well as two weak-IV confidence intervals, Anderson and Rubin (Anderson and Rubin 1949) confidence interval (Anderson and Rubin 1949) and the conditional likelihood ratio confidence interval (Moreira 2003). Finally, the code also conducts a sensitivity analysis if \eqn{Z_{i.}} is one-dimensional (i.e. there is only one instrument) using the method in Jiang et al. (2015).

Some procedures (e.g. conditional likelihood ratio test, sensitivity analysis with Anderson-Rubin) assume an additional linear model
\deqn{D_i = Z_{i.} \gamma +  \xi_i, E(\xi_i | X_{i.}, Z_{i.}) = 0}
}
\value{
  \code{ivmodel} returns an object of class "ivmodel".
  
  An object class "ivmodel" is a list containing the following components
  \item{alpha}{Significance level for the hypothesis tests}
  \item{beta0}{Null value of the hypothesis tests}
  \item{kClass}{A list from \code{kClassEst} function}
  \item{LIML}{A list from \code{LIMLEst} function}
  \item{Fuller}{A list from \code{FullerEst} function}
  \item{AR}{A list from \code{AR.test}}
  \item{CLR}{A list from \code{CLR.test}}
  In addition, if there is only one instrument, \code{ivreg} will generate an "ARsens" list within "ivmodel" object.
}
\references{
  Anderson, T. W. and Rubin, H. (1949). Estimation of the parameters of a single equation in a complete system of stochastic equations. \emph{Annals of Mathematical Statistics} 20, 46-63.\cr
  
  Fuller, W. (1977). Some properties of a modification of the limited information estimator. \emph{Econometrica}, 45, 939-953.  \cr
  
  Jiang, Y., Zhang, N and Small, D. (2015), Sensitivity analysis and power for instrumental variable studies, (under review of Biometrics).
  
  Moreira, M. J. (2003). A conditional likelihood ratio test for structural models. \emph{Econometrica} 71, 1027-1048.

  Sargan, J. D. (1958). The estimation of economic relationships using instrumental variables. Econometrica , 393-415. 
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\seealso{
See also \code{\link{kClassEst}}, \code{\link{LIMLEst}}, \code{\link{FullerEst}}, \code{\link{AR.test}}, and \code{\link{CLR.test}} for individual methods associated with \code{ivmodel}. For sensitivity analysis with the Anderson-Rubin test, see \code{\link{ARsens.test}}. \code{ivmodel} has \code{\link{print.ivmodel}}, \code{\link{confint.ivmodel}}, and \code{\link{coef.ivmodel}} methods associated with it.
}
\examples{
data(card.data)
# One instrument #
Y=card.data[,"lwage"]
D=card.data[,"educ"]
Z=card.data[,"nearc4"]
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661", "reg662", 
        "reg663", "reg664", "reg665", "reg666", "reg667", "reg668", "smsa66")
X=card.data[,Xname]
card.model1IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
card.model1IV

# Multiple instruments
Z = card.data[,c("nearc4","nearc2")]
card.model2IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
card.model2IV
}
\keyword{Instrumental variables}
\keyword{kClass estimation}
\keyword{Anderson-Rubin test}
\keyword{Sensitivity analysis with Anderson-Rubin test}
\keyword{Conditional likelihood ratio test}