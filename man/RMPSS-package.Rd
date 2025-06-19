\name{RMPSS-package}
\alias{RMPSS-package}
\alias{RMPSS}
\docType{package}
\title{
  Recursive Modified Pattern Search on Simplex (RMPSS)
}
\description{
 RMPSS is a black-box optimization tool which can be used to find the global minimum of
 any non-convex blackbox function where the parameter space belongs to an unit simplex, i.e.,
 we minimize f(x_1,...,x_n) where sum(x_i) = 1 and x_i>= 0 for all i.
}
\details{
  This package contains one main function RMPSS_opt. For more help, see '?RMPSS_opt' after
  loading the package into the local library.
}
\author{
Priyam Das.

Maintainer: Priyam Das <pdas@ncsu.edu>
}

\references{Das, Priyam,"Recursive Modified Pattern Search on High-Dimensional Simplex : A Blackbox Optimization Technique", Sankhya B 83 (Suppl 2), 440–483 (2021). https://doi.org/10.1007/s13571-020-00236-9  }

\keyword{Global optimization, Blackbox optimization, non-convex optimization, simplex constrained optimization. }
\seealso{

}
