#  File library/LDdiag/R/Copas.R
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


Copas <- function(obj, rate=0.5, ssize=1000) {
  scip <- options()$scipen
  options(scipen=4)
  c1 <- obj$call[1]
  c2 <- obj$call[2]
  c3 <- obj$call[3]
  c4 <- obj$call[4]
  calltype <- deparse(substitute(c1))
  calltype <- substr(calltype,1,nchar(calltype)-2)
  formula <- deparse(substitute(c2))
  formula <- substr(formula,1,nchar(formula)-2)
  if (calltype == "lm") {
    dt <- deparse(substitute(c3))
    dt <- substr(dt,1,nchar(dt)-2)
    data <- eval(substitute(dts,list(dts=as.name(dt))))
    Copas.lm(formula,data,rate,ssize)
    }
  else  {
    family <- c3 #deparse(substitute(c3))
    dt <- deparse(substitute(c4))
    dt <- substr(dt,1,nchar(dt)-2)
    dt <- eval(substitute(dts,list(dts=as.name(dt))))
    if(calltype=="glm") Copas.glm(formula, dt,family, rate,ssize)
    else cat("The function call doesn't support yet.\n\n")
  }
  options(scipen=scip)
}

Copas.lm <- function(formula, data, rate,ssize) {
  testf <- function(formula, data){
    fm <- formula
    sampos <- sample(1:nrow(data),round(nrow(data)*rate))
    dt1 <- data[sampos,]
    dt2 <- data[-sampos,]
    txt=paste("fit<- lm(", fm, ",", "dt1", ")", sep='')
    fit <- eval(parse(text=txt))
    yhat <- predict(fit, dt2,type="response")
    yname <- unlist(strsplit(formula, "~"))[1]
    yname <- substr(yname ,1, nchar(yname)-1)
    Y <- eval(substitute(dt2$vname,list(vname=as.name(yname))))
    Y <- Y - yhat
    newfit <- lm(Y ~ yhat)
    ses <- as.numeric(sqrt(diag(vcov(newfit))))
    tvals <- coef(newfit)/ses
    defs <-newfit$df
    pvalue <- 2*pt(abs(tvals),defs,lower.tail=FALSE)
    ab <- if (pvalue[1] <.05 && pvalue[2]<.05)  1
          else 0
    return(ab)
   }

  hold <- numeric(0)
  for (i in 1:ssize) hold[i] <- testf(formula, data)
  pt <- sum(hold)/ssize
  cat("\n  Of the",ssize,paste("times cross-validation (training set: ",rate*100,"% of the population). ",sep=""), paste(round(pt*100,1),"%",sep=""), "of the tests alpha!=0 and/or beta!=1.\n\n")
}



Copas.glm <- function(formula, data, family,rate, ssize) {
  testf <- function(formula, family, data){
    fm <- formula #as.name(deparse(substitute(formula)))
    sampos <- sample(1:nrow(data),round(nrow(data)*rate))
    dt1 <- data[sampos,]
    dt2 <- data[-sampos,]
    txt=paste("fit<- glm(", fm, ",",  "dt1,family=",  family, ")", sep='')
    fit <- eval(parse(text=txt))
    yhat <- predict(fit, dt2,type="response")
    yname <- unlist(strsplit(formula, "~"))[1]
    yname <- substr(yname ,1, nchar(yname)-1)
    Y <- eval(substitute(dt2$vname,list(vname=as.name(yname))))
    Y <- Y - yhat
    newfit <- lm(Y ~ yhat)
    ses <- as.numeric(sqrt(diag(vcov(newfit))))
    tvals <- coef(newfit)/ses
    defs <-newfit$df
    pvalue <- 2*pt(abs(tvals),defs,lower.tail=FALSE)
    ab <- if (pvalue[1] <.05 && pvalue[2]<.05)  1
          else 0
    return(ab)
  }
   hold <- numeric(0)
  for (i in 1:ssize) hold[i] <- testf(formula, family, data)
  pt <- sum(hold)/ssize
  cat("\n  Of the",ssize,paste("times cross-validation (training set: ",rate*100,"% of the population). ",sep=""), paste(round(pt*100,1),"%",sep=""), "of the tests alpha!=0 and/or beta!=1.\n\n")
}
