\name{linkpick}
\alias{linkpick}

\title{
Link Function Selection by Box-Cox Transformation
}

\description{
\code{linkpick} function selects candidate link functions by employing Box-Cox transformation. 
This function relies on \code{boxcox} funtion in \bold{MASS}. 
The dependend variable can only be positive numbers. If zeros exist, they should be replaced by small positive numbers before calling this function.\cr

The criteria used in this function are ("y" is assumed as dependent variable)
\itemize{
    \item[1.] if \eqn{\lambda = -2}, "inverse quadratic" link (\eqn{1/y^2}).
    \item[2.] if \eqn{\lambda = -1}, "inverse" link (1/y).
    \item[3.] if \eqn{\lambda = 0}, "log" link (log(y)).
    \item[4.] if \eqn{\lambda = 0.5}, "square root" link (\eqn{\sqrt{y} }).
    \item[5.] if \eqn{\lambda = 1}, "identity" (no transformation needed).
    \item[6.] if \eqn{\lambda = 2}, "square" link (\eqn{y^2}).
}
}

\usage{linkpick(formula,data)}

\arguments{
    \item{formula}{an formula as in \code{lm} function}
    \item{data}{data frame used for the model}
      }

\references{
Box, G.E.P. & Cox, D.R. (1964). An analysis of transformations. \emph{Journal of the Royal Statistical Society: Series B}. \bold{26}, 211-252

Venables, W. N. & Ripley, B.D. (2002). \emph{Modern Applied Statistics with S} (4th Edition). New York: Springer.
}

\examples{

quine$Days <- ifelse(quine$Days==0, 1, quine$Days)  #quine from MASS

linkpick(Days ~ ., data = quine)  #quine is from MASS

}

