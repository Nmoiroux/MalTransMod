
# predict-kdr with sex differential selection (on mating success for males and on number of oviposition events in female) # 
# selection on male before mating and on female after mating
predict_kdr <- function(q=0.1, W_F=c(1,1,1), W_M=c(1,1,1), N.gen=100){
  
  # initial population (hardy-weinberg equilibrium)
  p<- 1 - q
  
  RRf <- q^2    # frenquency RR females 
  SSf <- p^2    # frequency SS females
  RSf <- 2*p*q  # frequency RS females
  
  RRm <- RRf		# frenquency RR males
  SSm <- SSf		# frequency SS males
  RSm <- RSf		# frequency RS males
  
  # the data frame that tracks genotype frequencies over generations
  pop <- data.frame(RRf=RRf, RSf=RSf, SSf=SSf, RRm=RRm, RSm=RSm, SSm=SSm)
  
  for(i in 2:N.gen){
    
    RRf <- pop[i-1,1]		#frenquency RR females
    RSf <- pop[i-1,2]		#frequency RS females
    SSf <- pop[i-1,3]		#frequency SS females
    
    RRm <- pop[i-1,4] * W_M[1]		#frenquency RR males (taking into account relative mating succes)
    RSm <- pop[i-1,5] * W_M[2]		#frequency RS males 
    SSm <- pop[i-1,6] * W_M[3]		#frequency SS males 
    
    
    # next generation genotype probabilities (per female genotype)
    
    fRR_RR <- RRf * RRm +         # proportion of female RR carrying RR eggs (over total nb of females)
              0.5 * (RRf*RSm)
    fRR_RS <- RRf * SSm +         # female RR carrying RS eggs...
              0.5 * (RRf*RSm)
    
    fRS_RR <- 0.5 * (RRm*RSf) + 
              0.25* (RSm*RSf)
    fRS_RS <- 0.5 * (RRm*RSf) + 
              0.5 * (RSm*RSf) + 
              0.5 * (SSm*RSf)
    fRS_SS <- 0.25* (RSm*RSf) + 
              0.5 * (SSm*RSf)
    
    fSS_RS <- RRm * SSf + 
              0.5 * (SSf*RSm)
    fSS_SS <- SSf * SSm + 
              0.5 * (SSf*RSm)
    
    # next generation genotype quantities (adjusted for relative fitness of female = number of oviposition events)
    
    nextRR <- fRR_RR * W_F[1] + fRS_RR * W_F[2]
    nextRS <- fRR_RS * W_F[1] + fRS_RS * W_F[2] + fSS_RS * W_F[3]
    nextSS <- fRS_SS * W_F[2] + fSS_SS * W_F[3]
    
    # next generation normalised genotype frequencies
    
    pRR <- nextRR / sum(nextRR, nextRS, nextSS)
    pRS <- nextRS / sum(nextRR, nextRS, nextSS)
    pSS <- nextSS / sum(nextRR, nextRS, nextSS)
    
    # feed pop data frame
    pop[i,] <- c(pRR, pRS, pSS, pRR, pRS, pSS)
  } 
  
  # calculate kdr allelic frequency for all generations
  fkdr <- (2*pop$RRf + pop$RSf) / (2* (pop$RRf + pop$RSf + pop$SSf))
  return(fkdr) 
}

predict_kdr(W_F=c(1,0.8,0.66), W_M=c(0.67,1,0.55), N.gen=100)

# predict-kdr with sex differential selection (on mating success for males and on number of oviposition events in female) # tractable version
# selection on male and female occurs before mating
#' Title
#'
#' @param q 
#' @param W_F 
#' @param W_M 
#' @param N.gen 
#' @param sex.ratio
#'
#' @return
#' @export
#'
#' @examples
predict_kdr4 <- function(q=0.1, W_F=c(1,1,1), W_M=c(0.67,1,0.55), N.gen=100, sex.ratio = 0.5){
  
  # initialise first population (hardy-weinberg equilibrium)
  p <- 1 - q
  RRf <- q^2    # frenquency RR females 
  SSf <- p^2    # frequency SS females
  RSf <- 2*p*q  # frequency RS females
  
  RRm <- RRf		# frenquency RR males
  SSm <- SSf		# frequency SS males
  RSm <- RSf		# frequency RS males
  
  #
  pop <- data.frame(RRf=RRf, RSf=RSf, SSf=SSf, RRm=RRm, RSm=RSm, SSm=SSm)
  
  for(i in 2:N.gen){
    
    # genotype frequency after selection (generation i)
    RRf <- pop[i-1,1] * W_F[1]			#frenquency RR females
    RSf <- pop[i-1,2] * W_F[2]			#frequency RS females
    SSf <- pop[i-1,3] * W_F[3]			#frequency SS females
    
    RRm <- pop[i-1,4] * W_M[1]		  #frenquency RR males 
    RSm <- pop[i-1,5] * W_M[2]		  #frequency RS males 
    SSm <- pop[i-1,6] * W_M[3]		  #frequency SS males 
    
    # next generation genotype frequencies(random mating)
    mateRR <- RRf * RRm + 
      0.5 * (RRm*RSf) + 
      0.5 * (RRf*RSm) + 
      0.25* (RSm*RSf)
    
    mateRS <- RRf * SSm + 
      RRm * SSf + 
      0.5 * (RRm*RSf) + 
      0.5 * (RRf*RSm) + 
      0.5 * (RSm*RSf) + 
      0.5 * (SSm*RSf) + 
      0.5 * (SSf*RSm) 
    
    mateSS <- SSf * SSm  + 
      0.25* (RSm*RSf)  + 
      0.5 * (SSm*RSf)  + 
      0.5 * (SSf*RSm) 
    
    # next generation genotype probabilities (taking into account sex ratio)
    
    mRRf <- mRRm <- sex.ratio * mateRR
    mSSf <- mSSm <- sex.ratio * mateSS
    mRSf <- mRSm <- sex.ratio * mateRS
    
    # next generation normalised genotype frequencies per sex
    
    RRf2 <- mRRf / (mRRf + mSSf + mRSf)	#frenquency RR females
    SSf2 <- mSSf / (mRRf + mSSf + mRSf)	#frequency SS females
    RSf2 <- mRSf / (mRRf + mSSf + mRSf)	#frequency RS females
    
    RRm2 <- mRRm / (mRRm + mSSm + mRSm)	#frenquency RR males
    SSm2 <- mSSm / (mRRm + mSSm + mRSm)	#frequency SS males
    RSm2 <- mRSm / (mRRm + mSSm + mRSm)	#frequency RS males
    
    pop[i,] <- c(RRf2, RSf2, SSf2, RRm2, RSm2, SSm2)
  } 
  
  
  fkdr <- (pop$RRf + 1/2*pop$RSf) / (pop$RRf + pop$RSf + pop$SSf)
  return(fkdr) 
}

predict_kdr4(W_F=c(1,0.8,0.66), W_M=c(0.67,1,0.55), N.gen=100)


