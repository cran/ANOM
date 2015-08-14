\name{waterfilter}
\alias{waterfilter}
\docType{data}
\title{
Comparison of Water Filters
}
\description{
Filtering performances of seven brands of water filters, measured as the number of bacterial colonies growing on each device.
}
\usage{data(waterfilter)}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{brand}}{A factor with 7 levels giving the brands of water filters.}
    \item{\code{colonies}}{A numeric vector giving the number of bacterial colonies found on each filter.}
  }
}
\details{
A high number of bacterial colonies on a filter corresponds to good performance of this particular device. Note that the dataset is unbalanced (n=2 for brands 4 and 7, n=3 for all other brands).
}
\source{
Hsu, J. C. (1984) Ranking and selection and multiple comparisons with the best. In: Santner, T. J. and Tamhane, A. C. (Editors) Design of Experiments: Ranking and Selection (Essays in Honor of Robert E. Bechhofer). Marcel Dekker, New York, NY, pp. 23--33.
}
\references{
Westfall, P. H., Tobias, R. D., Wolfinger, R. D. (2011) Multiple Comparisons and Multiple Tests Using SAS, Second Edition. SAS Institute Inc., Cary, NC, pp. 592--593.
}
\examples{
data(waterfilter)
str(waterfilter)
}
\keyword{datasets}