\name{ANOM}
\alias{ANOM}
\title{
ANOM Decision Charts
}
\description{
Compute an analysis of means (i.e., a multiple contrast test involving comparisons of each group versus the grand mean) and draw a decision chart as commonly used in technometrics.
}
\usage{
ANOM(mc, xlabel=NULL, ylabel=NULL, printn=T, printp=T,
     stdep=NULL, stind=NULL, pst=NULL, pbin=NULL, bg="white")
}
\arguments{
  \item{mc}{An object of class \code{glht}, \code{SimCi}, \code{mctp}, or \code{binomRDci} involving group comparisons with the grand mean of all groups. See details.}
  \item{xlabel}{An optional character string specifying the label of the horizontal axis.}
  \item{ylabel}{An optional character string specifying the label of the vertical axis.}
  \item{printn}{A logical. Should per-group sample sizes be included in the chart? Default is \code{TRUE}.}
  \item{printp}{A logical. Should simultaneous p-values be included in the chart? Default is \code{TRUE}.}
  \item{stdep}{A numerical vector giving the values of the dependent variable. Only required if an object of class \code{SimCi} is inserted for \code{mc}, otherwise ignored. Default is \code{NULL}.}
  \item{stind}{A factor specifying the levels of the independent variable. Only required if an object of class \code{SimCi} is inserted for \code{mc}, otherwise ignored. Default is \code{NULL}.}
  \item{pst}{An object of class \code{SimTest}. Only required if an object of class \code{SimCi} is inserted for \code{mc} and simultaneous p-values are to be printed (\code{printp=TRUE}), otherwise ignored. Default is \code{NULL}.}
  \item{pbin}{An object of class \code{binomRDtest}. Only required if an object of class \code{binomRDci} is inserted for \code{mc} and simultaneous p-values are to be printed (\code{printp=TRUE}), otherwise ignored. Default is \code{NULL}.}
  \item{bg}{A character string. Should the plot's background be \code{"white"} (default) or \code{"gray"} (or \code{"grey"})?}
}
\details{
The 'standard' version of ANOM is invoked by inserting a \code{glht} object (created with function \code{glht} from package \code{multcomp} using a contrast matrix of type \code{GrandMean}) for \code{mc}. The \code{glht} object must be based on one of the model types \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, or \code{lmer}. That is, ANOM is feasible not only for simple linear (ANOVA) models with Gaussian data, but indeed for a broad range of (semi-)parametric models e.g., with Poisson or binomial data, hierarchical and clustered structures, and more (Hothorn et al. 2008). If the underlying model is a Poisson or binomial GLM (involving log and logit link functions, respectively), the effects are automatically transformed back to counts (Poisson) or proportions (binomial).

When analyzing binomial data, we need to distinguish between two data formats. As long as replicated observations of proportions are available for each group (i.e., a within-group variance can be estimated), we can fit a binomial GLM. However, if there is only one proportion per group (i.e., we have a 2-by-k data table), we need to fall back on some simpler procedure than a GLM. In the case of ANOM, we insert a \code{binomRDci} object (created with function \code{binomRDci} from package \code{MCPAN}) for \code{mc} with the contrasts being of type \code{GrandMean}. If printing out simultaneous p-values with the ANOM decision chart is desired (\code{printp=TRUE}), a \code{binomRDtest} object (generated with function \code{binomRDtest}) must be inserted for \code{pbin} (see examples). Notice that \code{printn=TRUE} prints the total sample size per group when inserting an object of class \code{binomRDci} for \code{mc}; by contrast, inserting a \code{glht} object based on a binomial GLM makes \code{printn=TRUE} print the number of independent observations of proportions per group.

There are two options for coping with heterogeneous variances. The 'standard' version using \code{glht} objects is basically designed for homoscedastic data but can also cope with unequal variances by employing sandwich covariance estimates (Herberich et al. 2010). All you have to do is set the option \code{vcov=vcovHC} in the \code{glht} call (see examples).

A different heterogeneity adjustment using multiple Satterthwaite degrees of freedom (Hasler and Hothorn 2008) is invoked by inserting a \code{SimCi} object (created with function \code{SimCiDiff} from package \code{SimComp} with the setting \code{covar.equal=FALSE}) for \code{mc}. Notice that \code{covar.equal=TRUE} would compute a 'standard' multiple contrast test assuming equal variances, which is just the opposite of what we want.

Another instance when the \code{SimComp} package proves useful is for ANOM with ratios i.e., we assess each group's percentage change in comparison to the grand mean (which is always 100\%). Again, we insert an object of class \code{SimCi} for \code{mc} but now we generated it with function \code{SimCiRatio} instead of \code{SimCiDiff}. Setting \code{covar.equal} to \code{TRUE} gives the homoscedastic version whereas \code{FALSE} calls the heteroscedastic variant using multiple degrees of freedom. 

All ANOM functionality based on objects from \code{SimComp} requires you to submit the data (dependent and independent variable via \code{stdep} and \code{stind}, respectively) separately (see examples). Moreover, if you wish to print simultaneous p-values with your ANOM decision chart, you may compute them using function \code{SimTestDiff} or \code{SimTestRat} and insert the resulting object for \code{pst}. Again, don't forget to set the option \code{covar.equal=FALSE} if you want to account for heteroscedasticity. Make sure in all cases to perform comparisons of type \code{GrandMean}.

Nonparametric ANOM is performed by inserting an object of class \code{mctp} (created with function \code{mctp} from package \code{nparcomp}) for \code{mc}. It is based on estimation of relative effects via global pseudo-rankings as proposed by Konietschke et al. (2012). This procedure can naturally cope with heteroscedasticity in the data.

The relative effect of two independent random variables X_1 and X_2 following some distributions F_1 and F_2, respectively, is generally defined as
\deqn{p = P(X_1<X_2) + 0.5 P(X_1=X_2).}
Loosely speaking, p is the probability that X_1 takes smaller values than X_2 (plus half the probability of taking equal values). Hence when p<0.5, X_1 is stochastically more likely to take larger values than X_2, and vice versa for p>0.5.

Creating the \code{mctp} object requires to specify a grand-mean-type contrast matrix by hand (see examples) as it is not among the options provided by the package \code{nparcomp}. For the test statistics you may choose one out of three asymptotic approximation methods:
\itemize{
\item a multivariate t approximation with Satterthwaite degrees of freedom (\code{asy.method="mult.t"}), which works reasonably well most of the time,
\item a multivariate normal approximation (\code{asy.method="normal"}), which is unfavorable with small sample sizes,
\item a Fisher transform (\code{asy.method="fisher"}), which ensures that the decision limits preserve the range of [-1, 1].
}
Make sure that the argument \code{correlation} in the function \code{mctp} is set to \code{TRUE}. Do not use the function \code{nparcomp} instead of \code{mctp} since it does not involve global ranking and is inoperative with ANOM-type contrast matrices.
}
\value{
An ANOM decision chart.
}
\references{
Djira, G. D., Hothorn, L. A. (2009) Detecting relative changes in multiple comparisons with an overall mean. Journal of Quality Technology 41(1), 60-65.

Hasler, M. and Hothorn, L. A. (2008) Multiple contrast tests in the presence of heteroscedasticity. Biometrical Journal 50(5), 793--800.

Herberich, E., Sikorski, J., Hothorn, T. (2010) A robust procedure for comparing multiple means under heteroscedasticity in unbalanced designs. PLoS One 5(3), e9788.

Hothorn, T., Bretz, F., Westfall, P. (2008) Simultaneous inference in general parametric models. Biometrical Journal 50(3), 346--363.

Konietschke, F., Hothorn, L. A., Brunner, F. (2012) Rank-based multiple test procedures and simultaneous confidence intervals. Electronic Journal of Statistics 6, 738--759.

Pallmann, P. and Hothorn, L. A. (2016) Analysis of means (ANOM): A generalized approach using R. Journal of Applied Statistics. DOI: 10.1080/02664763.2015.1117584.
}
\author{
Philip Pallmann \email{p.pallmann@lancaster.ac.uk}
}
\note{
Notice that some of the more sophisticated ANOM variants (ratios, nonparametric, heteroscedastic with multiple degrees of freedom) are limited to one-way layouts without covariates.
}
\examples{\donttest{
###############################################
### Standard ANOM (Gaussian, homoscedastic) ###
###############################################

### Devices of which brand filter bacteria significantly worse?
head(waterfilter)
str(waterfilter)

library(multcomp)
model <- lm(colonies ~ brand, waterfilter)
hom <- glht(model, mcp(brand="GrandMean"), alternative="less")
ANOM(hom)

############################
### Heteroscedastic ANOM ###
############################

## With sandwich covariance matrix estimate (Herberich et al. 2010)

library(multcomp)
library(sandwich)
het1 <- glht(model, mcp(brand="GrandMean"), alternative="less", vcov=vcovHC)
ANOM(het1)

## With multiple degrees of freedom (Hasler and Hothorn 2008)

library(SimComp)
het2 <- SimCiDiff(data=waterfilter, grp="brand", resp="colonies",
                  type="GrandMean", alternative="less", covar.equal=FALSE)
het2p <- SimTestDiff(data=waterfilter, grp="brand", resp="colonies",
                     type="GrandMean", alternative="less", covar.equal=FALSE)
ANOM(het2, stdep=waterfilter$colonies, stind=waterfilter$brand, pst=het2p)

#######################
### ANOM for ratios ###
#######################

## Homoscedastic

library(SimComp)
rel <- SimCiRat(data=waterfilter, grp="brand", resp="colonies",
                type="GrandMean", alternative="less", covar.equal=TRUE)
relp <- SimTestRat(data=waterfilter, grp="brand", resp="colonies",
                   type="GrandMean", alternative="less", covar.equal=TRUE)
ANOM(rel, stdep=waterfilter$colonies, stind=waterfilter$brand, pst=relp)

## Heteroscedastic (with multiple degrees of freedom)

library(SimComp)
relh <- SimCiRat(data=waterfilter, grp="brand", resp="colonies",
                 type="GrandMean", alternative="less", covar.equal=FALSE)
relhp <- SimTestRat(data=waterfilter, grp="brand", resp="colonies",
                    type="GrandMean", alternative="less", covar.equal=FALSE)
ANOM(relh, stdep=waterfilter$colonies, stind=waterfilter$brand, pst=relhp)

##########################
### Nonparametric ANOM ###
##########################

# Compute sample sizes per group
ss <- tapply(waterfilter$colonies, waterfilter$brand, length)
# Build a grand-mean-type contrast matrix
library(multcomp)
Mat <- contrMat(ss, "GrandMean")

## Using a multivariate t approximation

library(nparcomp)
mult <- mctp(colonies ~ brand, data=waterfilter, type="UserDefined",
             contrast.matrix=Mat, alternative="less", info=FALSE,
             correlation=TRUE, asy.method="mult.t")
ANOM(mult)
           
## Using a range-preserving Fisher transform

library(nparcomp)
fish <- mctp(colonies ~ brand, data=waterfilter, type="UserDefined",
             contrast.matrix=Mat, alternative="less", info=FALSE,
             correlation=TRUE, asy.method="fisher")
ANOM(fish)

#####################################
### ANOM for binomial proportions ###
#####################################

### Which schools' math achievements differ from the grand mean?

head(math)
str(math)

## Based on Wald-type confidence intervals

library(MCPAN)
wald <- binomRDci(n=math$enrolled, x=math$proficient, names=math$school,
                  alternative="two.sided", method="Wald", type="GrandMean")
waldp <- binomRDtest(n=math$enrolled, x=math$proficient, names=math$school,
                     alternative="two.sided", method="Wald", type="GrandMean")
ANOM(wald, pbin=waldp)

## Based on add-2 confidence intervals

library(MCPAN)
add2 <- binomRDci(n=math$enrolled, x=math$proficient, names=math$school,
                  alternative="two.sided", method="ADD2", type="GrandMean")
add2p <- binomRDtest(n=math$enrolled, x=math$proficient, names=math$school,
                     alternative="two.sided", method="ADD2", type="GrandMean")
ANOM(add2, pbin=add2p)

##########################
### ANOM for variances ###
##########################

### Springs of which brand are significantly more variable?

head(spring)
str(spring)

# Compute the median weight per brand
spring$median <- tapply(spring$weight, spring$brand, median)[spring$brand]
# Compute the absolute deviations from the median (robust Levene residuals)
spring$absdev <- with(spring, abs(weight - median))

library(multcomp)
mod <- lm(absdev ~ brand, spring)
test <- glht(mod, mcp(brand="GrandMean"))
ANOM(test)
}}
\keyword{~dplot}