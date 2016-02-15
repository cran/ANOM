\name{ANOMgen}
\alias{ANOMgen}
\title{
Generic Function for Drawing ANOM Decision Charts
}
\description{
Graphical representation of the analysis of means: convert simultaneous confidence intervals (that were computed with ANY method) into ANOM decision limits and draw a decision chart as commonly used in technometrics.
}
\usage{
ANOMgen(mu, n=NULL, gm=NULL, lo, up, names, alternative="two.sided",
        xlabel="Group", ylabel="Endpoint", printn=T, p=NULL, bg="white")
}
\arguments{
  \item{mu}{A numeric vector of group means.}
  \item{n}{A numeric vector of sample sizes per group. Either \code{n} or \code{gm} must be provided.}
  \item{gm}{A single numeric value giving the grand mean of all groups. Either \code{n} or \code{gm} must be provided.}
  \item{lo}{A numeric vector of lower (simultaneous) confidence interval bounds for comparisons to the grand mean.}
  \item{up}{A numeric vector of upper (simultaneous) confidence interval bounds for comparisons to the grand mean.}
  \item{names}{An (optional) vector of characters specifying the groups' names.}
  \item{alternative}{A character string indicating the direction of the alternative hypothesis. Default is \code{"two.sided"}, but may be changed to one-sided alternatives (either \code{"greater"} or \code{"less"}).}
  \item{xlabel}{A character string specifying the label of the horizontal axis.}
  \item{ylabel}{A character string specifying the label of the vertical axis.}
  \item{printn}{A logical. Should per-group sample sizes be included in the chart? Default is \code{TRUE}. If \code{n} ist left at \code{NULL}, the function automatically sets \code{printn} to \code{FALSE}.}
  \item{p}{An (optional) numeric vector of (simultaneous) p-values to be printed.}
  \item{bg}{A character string. Should the plot's background be \code{"white"} (default) or \code{"gray"} (or \code{"grey"})?}
}
\details{
This is a generic tool that translates (simultaneous) confidence intervals into ANOM decision limits.
}
\value{
An ANOM decision chart.
}
\references{
Pallmann, P. and Hothorn, L. A. (2016) Analysis of means (ANOM): A generalized approach using R. Journal of Applied Statistics. DOI: 10.1080/02664763.2015.1117584.
}
\author{
Philip Pallmann \email{p.pallmann@lancaster.ac.uk}
}
\note{
The confidence intervals must arise from comparisons to the grand mean; otherwise the ANOM chart is meaningless!
}
\examples{\donttest{
### A toy example (n given, two-sided)
groupmeans <- c(2.8, 2.3, 3.4, 5.6)
samplesizes <- c(5, 5, 10, 5)
low <- c(-1.2, -1.7, -0.4, 1.6)
upp <- c(-0.2, -0.7, 0.2, 2.6)
names <- c("1st", "2nd", "3rd", "4th")
ANOMgen(mu=groupmeans, n=samplesizes, lo=low, up=upp, names=names, alternative="two.sided")

### Another toy example (gm given, one-sided, with p-values)
groupmeans <- c(2.8, 2.3, 3.4, 5.6)
gm <- 3.5
low <- rep(-Inf, 4)
upp <- c(-0.2, -0.7, 0.2, 2.6)
names <- c("1st", "2nd", "3rd", "4th")
pvalues <- c(0.01, 0.003, 0.8, 1)
ANOMgen(mu=groupmeans, gm=gm, lo=low, up=upp, names=names, alternative="less", p=pvalues)
}}
\keyword{~dplot}