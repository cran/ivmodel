\name{ivmodel-internal}
\title{Internal \code{ivmodel} functions}
\alias{invTwobyTwoSymMatrix}
\alias{condPvalue}
\alias{qrrank}
\alias{qrRM}
\alias{quadSolver}
\description{Internal \code{ivmodel} functions}
\usage{
  invTwobyTwoSymMatrix(M)
  quadSolver(a,b,c)
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\details{
  \code{invTwobyTwoSymMatrix} is a helper function to invert two-by-two symmetric matrices. \code{condPvalue} is a helper function for \code{CLR.test}.
  \code{qrrank} and \code{qrRM} are helper functions for the sparse QR package.
  \code{quadSolver} is a helper function to invert to solve quadratic equations.
}
\keyword{internal}