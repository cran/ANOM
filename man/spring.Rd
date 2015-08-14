\name{spring}
\alias{spring}
\docType{data}
\title{
Stiffness of Springs
}
\description{
Weights required to stretch springs of four brands by 0.1 inches.
}
\usage{data(spring)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{\code{brand}}{A factor with 4 levels giving the brands of springs.}
    \item{\code{weight}}{A numeric vector giving the weight required to extend the spring by 0.1 inches.}
  }
}
\source{
Nelson, P. R., Wludyka, P. S., Copeland, K. A. F. (2005) The Analysis of Means: A Graphical Method for Comparing Means, Rates, and Proportions. Society for Industrial and Applied Mathematics (SIAM), Philadelphia, PA, and American Statistical Association (ASA), Alexandria, VA, p. 53.
}
\examples{
data(spring)
str(spring)
}
\keyword{datasets}
