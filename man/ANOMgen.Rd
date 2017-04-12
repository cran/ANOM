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
        xlabel="Group", ylabel="Endpoint", printn=TRUE, p=NULL, bg="white",
        bgrid=TRUE, axlsize=18, axtsize=25, npsize=5, psize=5, lwidth=1,
        dlstyle="dashed", fillcol="darkgray")
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
  \item{bgrid}{A logical. Should background grid lines be plotted? Default is \code{TRUE}.}
  \item{axlsize}{A numerical value specifying the font size of the axis labels. Default is \code{18}.}
  \item{axtsize}{A numerical value specifying the font size of the axis titles. Default is \code{25}.}
  \item{npsize}{A numerical value specifying the font size of the sample sizes and p-values (if printed). Default is \code{5}.}
  \item{psize}{A numerical value specifying the size of the points (group means). Default is \code{5}.}
  \item{lwidth}{A numerical value specifying the width of the lines (grand mean, decision limits, vertical connections). Default is \code{1}.}
  \item{dlstyle}{A character string specifying the style of the decision limits. Default is \code{"dashed"}.}
  \item{fillcol}{A character string specifying the color of the area of no significant deviation from the grand mean. Default is \code{"darkgray"}.}
}
\details{
This is a generic tool that translates (simultaneous) confidence intervals into ANOM decision limits.
}
\value{
An ANOM decision chart.
}
\references{
Pallmann, P. and Hothorn, L. A. (2016) Analysis of means (ANOM): A generalized approach using R. Journal of Applied Statistics, 43(8), 1541--1560.
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