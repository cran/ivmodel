\name{AR.size}
\alias{AR.size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{AR.size
%%  ~~function to do ... ~~
}
\description{\code{AR.size} computes the minimum sample size required for achieving certain power of Anderson-Rubin (1949) test.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AR.size(power, k, l, beta, gamma, Zadj_sq, sigmau, sigmav, rho, 
        alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{The goal of power achieving over a constant.
%%     ~~Describe \code{power} here~~
}
  \item{k}{Number of exogenous variables.
%%     ~~Describe \code{k} here~~
}
  \item{l}{Number of instrumental variables.
%%     ~~Describe \code{l} here~~
}
  \item{beta}{True causal effect minus null hypothesis causal effect.
%%     ~~Describe \code{beta} here~~
}
  \item{gamma}{Regression coefficient for effect of instruments on treatment.
%%     ~~Describe \code{gamma} here~~
}
  \item{Zadj_sq}{Variance of instruments after regressed on the observed variables.
%%     ~~Describe \code{Zadj_sq} here~~
}
  \item{sigmau}{Standard deviation of potential outcome under control. (structural error for y)
%%     ~~Describe \code{sigmau} here~~
}
  \item{sigmav}{Standard deviation of error from regressing treatment on instruments.
%%     ~~Describe \code{sigmav} here~~
}
  \item{rho}{Correlation between u(potential outcome under control ) and v (error from regressing treatment on instrument).
%%     ~~Describe \code{rho} here~~
}
  \item{alpha}{Significance level.
%%     ~~Describe \code{alpha} here~~
}
}

\value{Minimum sample size required for achieving certain power of Anderson-Rubin (1949) test.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.
%% ~put references to the literature/web site here ~
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
### Minimum sample size needed for power of sensitivity analysis over 0.8. 
### In a study where the null hypothesis causal effect is 0,
### the true causal effect is 1, the sample size is 250, the instrument is
### binary with probability .5 (so variance = .25), the standard deviation
### of potential outcome under control is 1, the effect of the instrument
### is to increase the probability of a binary treatment being 1 from .25 to 
### .75.  The function sigmav.func computes the SD of v for a binary insrument, 
### binary treatment. The correlation between u and v is assumed to be .5.  The 
### significance level for the study will be alpha = .05
### the calculated sigmav is .4330127

### minimum sample size for Anderson-Rubin test
AR.size(power=0.8, k=1, l=1, beta=1, gamma=.5, Zadj_sq=.25, sigmau=1, 
        sigmav=.4330127, rho=.5, alpha = 0.05)

}
