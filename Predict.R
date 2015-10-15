library(survival)
library("survcomp")

Predict <- function(features.train, features.test, Surv.train, Surv.test, Cens.train, Cens.test)
{
  out.train <- c()
  out.test <- c()
  for(i in 1 : ncol(features.train))
  {
    if (length(which(table(features.train[,i])> 0.8*length(Surv.train)))>0)  # discard the flat values.
    {
      out.train <- c(out.train,i)    
    }
    if (length(which(table(features.test[,i])> 0.8*length(Surv.test)))>0)    # discard the flat values.
    {
      out.test <- c(out.test,i)    
    }
  }
  if (length(out.train)>0 | length(out.test)>0)
  {
    out <- union(out.train,out.test)
    features.train <- features.train[,-out]
    features.test <- features.test[,-out]
  }
  
  model <- coxph(Surv(Surv.train,Cens.train,type='right') ~ features.train, control = coxph.control(iter.max = 500))
  model$coefficients[is.na(model$coefficients)]=0           # convert NAs to zeros if any.
  model.predict <- as.matrix(features.test) %*% as.numeric(model$coefficients)
  c.index <- concordance.index(model.predict, Surv.test, Cens.test, method="noether")
  return(c.index)
}
