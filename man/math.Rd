\name{math}
\alias{math}
\docType{data}
\title{
Math Proficiency Scores
}
\description{
Proportion of fifth graders with proficient math test scores in 10 elementary schools.
}
\usage{data(math)}
\format{
  A data frame with 10 observations on the following 3 variables.
  \describe{
    \item{\code{school}}{A factor with 10 levels giving the ID of the school.}
    \item{\code{enrolled}}{A numeric vector giving the number of students taking part in the math test.}
    \item{\code{proficient}}{A numeric vector giving the number of students with proficient math test scores.}
  }
}
\details{
A study compared math achievements of students from 10 elementary schools in a U.S. district; 6 of them were conventional neighborhood schools (\code{N1}--\code{N6}) and 4 alternative schools (\code{A1}--\code{A4}). 563 fifth graders took standardized math tests, and each school's proportion of students who scored proficient was recorded.
}
\source{
Nelson, P. R., Wludyka, P. S., Copeland, K. A. F. (2005) The Analysis of Means: A Graphical Method for Comparing Means, Rates, and Proportions. Society for Industrial and Applied Mathematics (SIAM), Philadelphia, PA, and American Statistical Association (ASA), Alexandria, VA, pp. 42--43.
}
\examples{
data(math)
str(math)
}
\keyword{datasets}