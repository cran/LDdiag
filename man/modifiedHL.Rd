

\name{modifiedHL}
\alias{modifiedHL}

\title{Modified Hosmer-Lemeshow Diagnostic Test}

\description{
The Hosmer-Lemeshow Statistic is a measure of goodness of
fit. According to Hosmer and Lemeshow, this function 1)partitions the observations into 10 equal sized groups
according to their predicted values, 2) runs a F-test for equality of mean residuals of the 10 groups.
If null hypothesis is rejected, the model doesn't fit the data well.
}

\usage{modifiedHL(object,digits=4)}

\arguments{
    \item{object}{an object created by \code{lm} or \code{glm} function}
    \item{digits}{round to decimals. Default is 4}
      }

\references{
Hosmer, D.W., Lemeshow,  S.(2000).  \emph{Applied Logistic Regression} (2nd Edition).  New York: Wiley. \cr
}

\examples{

ex <- lm(Days ~ ., data = quine)  #quine is from MASS
modifiedHL(ex)

quine$Days <- ifelse(quine$Days==0, 1, quine$Days)
ex <- glm(Days ~ ., family = Gamma(link="log"), data = quine)
modifiedHL(ex)

}

