Score <- function(Cox.Uni.beta, Cox.Uni.p, Spearman.beta, Spearman.p, Genes.Sig, HR.o=0.5, q=0.05)
{

  Spearman.q <- matrix(1, nrow(Spearman.beta), ncol(Spearman.beta))
  for (i in 1 : ncol(Spearman.beta))
  {
    Spearman.q[,i] <- p.adjust(Spearman.p[,i], method="BH")
  }  
  q.in <- (Spearman.q <= q)
  sig.q.num <- colSums(q.in)  # Number of associated GE features per each variant feature
  zq <- which(sig.q.num > 0)  # which variant features have association with the GE features
  
  S.Obs <- c()        # observed score
  S.pvalue <- c()     # P-value of the significance of each score
  Ind.list <- list()  # indices of the selected GE features per each significant score S_n.The number of
                      # entries of Ind.list equals the number of variants with significant S_n as in Equation (2)
  j <- 1
  
  for (i in 1 : length(zq))
  {
    ind <- which(q.in[,zq[i]])                          # indices of the associated GE features per variant
    Cox.pvalue <- Cox.Uni.p[Genes.Sig[ind]]
    Cox.beta <- Cox.Uni.beta[Genes.Sig[ind]]
    HR <- exp(abs(Cox.beta))
    Rho <- Spearman.beta[ind,zq[i]]
    S.g <- vector(length=length(ind), mode="numeric")
    S.g <- (HR - HR.o) - (1 - abs(Rho))                 # Equation (1) in the paper
    S <- sum(S.g[S.g>0])    
    if (S > 0)
    {
      S <- S / sum(S.g>0)
      temp <- which(S.g>0)
      Ind.list <- c(Ind.list, list(Genes.Sig[ind[temp]]))
      # if all S.g's <= 0, Ind.list will have an empy entry
      # Ind.list has the original indices of mRNA features
      # Gene_Sig: significant mRNA features
      # ind: part of Gene.Sig that is correlated with the variant
      # temp: part of ind that has positive score
      
      S.Obs[j] <- S
      
      # Permutation test to construct S.Null and estimate P-value
      S.Null <- matrix(0, 1000, 1)
      for (k in 1 : 1000)
      {
        HR.ind <- sample(length(ind))    # Randome permutaion to the HR values
        Rho.ind <- sample(length(ind))   # Randome permutaion to the Rho values
        HR.temp <- HR[HR.ind]
        Rho.ind <- Rho[Rho.ind]
        
        S.g.temp <- (HR.temp - HR.o) - (1 - abs(Rho.ind))
        S.Null[k] <- sum(S.g.temp[S.g.temp>0])
        if(S.Null[k] > 0)
        {
          S.Null[k] <- S.Null[k] / sum(S.g.temp>0)
        }
      }
      S.pvalue[j] <- sum(S.Null >= S) / 1000
      j <- j + 1
    }    
  }
  
  S.fdr <- p.adjust(S.pvalue, method="BH")

  S.fdr.sig <- which(S.fdr <= 0.05)  # The significant S_n's
  
  Ind.list <- Ind.list[S.fdr.sig]

  return(Ind.list)
}
