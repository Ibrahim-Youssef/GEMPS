library(survival)
source('~/RStudio/Score.R')
source('~/RStudio/Predict.R')


Survival_Analysis_GEMPS <- function(
  GE,        # a matrix of M Patients x N mRNA features
  train.ind, # a vector specifying the indices of the samples used for training the Cox model
             # the testing indices will be generated automatically by the code from the remaining samples
  Surv,      # a vector of length M for the survival data
  Cens,      # a vector of length M for the censoring data
  Clinical,  # a matrix of M patients x L clinical features
  Variant,   # a matrix of M patients x V variant features
  HR.o=0.5,  # value of the HR.o parameter in Equation (1) in the paper
  q=0.05     # the FDR threshold. Default is 0.05
  )

  
{
  NPatients <- nrow(GE)       # number of the patients (samples)
  NClinical <- ncol(Clinical) # number of the clinical features
  Out <- c()
  for(i in 1 : ncol(GE))
  {
    if (length(which(table(GE[,i])> 0.8*NPatients))>0) # discard the flat values.
    {
      Out <- c(Out,i)    
    }
  }
  if (length(Out) > 0)
  {
    GE <- GE[,-Out]
  }
  GE.original.ind <- setdiff(c(1:ncol(GE)), Out)
  
  Out <- c()
  for(i in 1 : ncol(Variant))
  {
    if (length(which(table(Variant[,i])> 0.8*NPatients))>0) # discard the flat values.
    {
      Out <- c(Out,i)    
    }
  }
  if (length(Out) > 0)
  {
    Variant <- Variant[,-Out]
  }
  
  test.ind <- setdiff(c(1:NPatients), train.ind) # generating the indices of the testing samples
  
  GE.train <- GE[train.ind,]
  GE.test <- GE[test.ind,]
  
  Surv.train <- Surv[train.ind]
  Surv.test <- Surv[test.ind]
  Cens.train <- Cens[train.ind]
  Cens.test <- Cens[test.ind]
  
  Clinical.train <- Clinical[train.ind,]
  Clinical.test <- Clinical[test.ind,]
  
  Variant.train <- Variant[train.ind,]
  Variant.test <- Variant[test.ind,]
  
  # The univariate Cox model
  UCox_Gene_p <- matrix(1, ncol(GE.train), (NClinical+1))  # Covariates per time are a GE feature plus the clinical features
  UCox_Gene_beta <- matrix(0, ncol(GE.train), (NClinical+1))
  for(i in 1 : ncol(GE.train))
  {
    Model <- coxph(Surv(Surv.train,Cens.train,type='right') ~ GE.train[,i] + Clinical.train, control = coxph.control(iter.max = 500))  
    UCox_Gene_p[i,] <- summary(Model)$coefficients[,5] # Wald Test
    UCox_Gene_beta[i,] <- as.numeric(Model$coefficients)
  }
  GE_Sig <- which(UCox_Gene_p[,1] <= 0.05) # Indices of the significant mRNA features based on the univariate Cox model
  
  
  # The eQTL-like process between variants and the gene expression phenotypes using the Spearman's coefficient
  Spearman_p <- matrix(1, length(GE_Sig), ncol(Variant.train))
  Spearman_beta <- matrix(0, length(GE_Sig), ncol(Variant.train))
  for(i in 1 : length(GE_Sig))
  {
    for(j in 1 : ncol(Variant.train))
    {
      Correlation <- cor.test(GE.train[,GE_Sig[i]], Variant.train[,j], method = "spearman")
      Spearman_p[i,j] <- Correlation$p.value
      Spearman_beta[i,j] <- Correlation$estimate
    }    
  }  
  
  out.GEMPS <- c()
  # Computing the GEMPS Score in Equations (1) and (2) and checking for significance using the function "Score"
  L <- Score(UCox_Gene_beta, UCox_Gene_p, Spearman_beta, Spearman_p, GE_Sig, HR.o, q)
  if(length(L) > 0)
  {
    # The union operation to include all genes associated with individual variants or having polymorphic forms of variants
    genes.ind.GEMPS <- c()
    for (i in 1 : length(L))
    {
      genes.ind.GEMPS <- union(genes.ind.GEMPS, L[[i]])
    }
    
    genes.ind.GEMPS <- sort(genes.ind.GEMPS)
    GE.train.GEMPS <- as.matrix(GE.train[,genes.ind.GEMPS])
    GE.test.GEMPS <- as.matrix(GE.test[,genes.ind.GEMPS])
    cindex.GEMPS <- try(Predict(GE.train.GEMPS, GE.test.GEMPS, Surv.train, Surv.test, Cens.train, Cens.test)$c.index)
    genes.ind.GEMPS.original <- GE.original.ind[genes.ind.GEMPS]
    
    if(class(cindex.GEMPS) != "try-error")
    {
      out.GEMPS <- c(list(cindex.GEMPS), list(genes.ind.GEMPS.original))
    }      
  }
  return(out.GEMPS)
}