\name{hemoglobin}
\alias{hemoglobin}
\docType{data}
\title{
Hemoglobin Levels
}
\description{
Hemoglobin levels of 30 male cancer patients treated with radiation or chemotherapy and one of three drugs.
}
\usage{data(hemoglobin)}
\format{
  A data frame with 30 observations on the following 3 variables.
  \describe{
    \item{\code{therapy}}{A factor with 2 levels giving the types of therapy.}
    \item{\code{drug}}{A factor with 3 levels giving the drugs administered.}
    \item{\code{level}}{A numeric vector giving the patients' hemoglobin levels.}
  }
}
\details{
This is a complete balanced two-way layout. 15 male cancer patients were radiated, and another 15 underwent chemotherapy. In addition, the patients were treated with either drug 1, 2, or 3. The endpoint of interest was the level of hemoglobin (in grams per deciliter blood).
}
\source{
Nelson, P. R., Wludyka, P. S., Copeland, K. A. F. (2005) The Analysis of Means: A Graphical Method for Comparing Means, Rates, and Proportions. Society for Industrial and Applied Mathematics (SIAM), Philadelphia, PA, and American Statistical Association (ASA), Alexandria, VA, pp. 71 ff.
}
\examples{
data(hemoglobin)
str(hemoglobin)
}
\keyword{datasets}
