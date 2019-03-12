#' `VLAIB` return the average number of infectious bites that a vector gives according to a behavior and mortality model
#'
#' @param nsim unused but required argument, for compatibility / comparison with other function (ex: VLAIB_IBM)
#' @param S baseline survival rate
#' @param g duration of gonotrophic cycle (positive integer, default = 3)
#' @param Du Diversion probability when entering a hut without LLIN (default = 0.43)
#' @param Dp Diversion probability when entering a hut with LLIN (default = 0.3)
#' @param m1u pre-bite feeding related mortality probability when faced to un unprotected host (default = 0.03)
#' @param m1p pre-bite feeding related mortality probability when faced to un LLIN protected host (default = 0.72)
#' @param m2u post-bite mortality probability when feeding on an unprotected host (default = 0.005)
#' @param m2p post-bite mortality probability when feeding on an LLIN protected host (default = 0.21)
#' @param Nh Number of humans in the community (default = 1000)
#' @param Uh proportion of humans that use LLINs (default = 0.6)
#' @param pi proportion of exposure to bite that occurs during which LLIN is in use (default = 0.9)
#' @param Pllin preference for LLIN protected human (against uprotected human) as recorded in a dual choice olfactometer (default = 0.5, i.e inert LLIN)
#' @param k infectiousness: probability that a vector become infected (exposed) while taking a blood meal on an infectious host (default = 0.1)
#' @param n duration of extrinsic incubation period of the parasite (in days, positive integer, default = 11)
#' @param Ih Plasmodium falciparum prevalence rate in the human population (default = 0.5)
#'
#' @return a vector of several output elements including the VLAIB. Each element of the returning vector is named:
#' "Eu", "Ep", "Pd", "Sd", "Pf_u", "Pf_p", "Pf", "F_u", "Sf", "Sl", "PfA", "BA","OvA", "Pl", "VLAIB"  
#' @export
#'
#' @examples
#' VLAIB()
#' VLAIB(Ih = 0.1)
#' VLAIB()["VLAIB"]
VLAIB <- function(nsim=1000, S = 0.9, g = 3, Du = 0.43, Dp= 0.3,m1u = 0.05, m1p = 0.72, m2u = 0.005, m2p = 0.21, 
									Nh = 1000, Uh = 0.6, pi = 0.9, Pllin = 0.5, k = 0.1, n = 11, Ih = 0.5){
  
  ## proba of encountering an unprotected or LLIN protected human
  Np <- Uh*pi*Nh                        # average number of people protected by an LLIN - Expression (1)
  Nu <- Nh - Np                        # average number of unprotected people - Expression (2)
  Cpp <- dhyper(2, round(Np), round(Nu), 2, log = FALSE) # probability that a host-seeking Anopheles will be faced to a choice between two LLIN protected hosts - Expression (4)
  Cpu <- dhyper(1, round(Np), round(Nu), 2, log = FALSE) # ... choice between a LLIN protected host and an unprotected host - Expression (5)
  Cuu <- dhyper(0, round(Np), round(Nu), 2, log = FALSE) # ... choice between two unprotected hosts- Expression (6)
  Eu <- Cuu + Cpu * (1 - Pllin)				  # proba to enter a house with an unprotected human - Expression (7)
  Ep <- Cpp + Cpu * Pllin			          # proba to enter a house with a protected human - Expression (8)
  
  
  ## diversion, feeding and mortality probabilities
	fi1_u <- 1 - m1u 		      # Successful feeding probability when ettempting to feed (non diverted mosquitoes) in a hut without LLIN - Expression (9)
  fi1_p <- 1 - m1p          # Successful feeding probability when ettempting to feed (non diverted mosquitoes) in a hut with LLIN- Expression (12)
  fi_u  <- (1-Du) * fi1_u   # Successful feeding probability when entering a hut with an unprotected human - Expression (10)
  fi_p  <- (1-Dp) * fi1_p   # Successful feeding probability when entering a hut with an LLIN protected human- Expression (13)
  S2u <- 1 - m2u            # post-bite survival in hut with unprotected people - Expression (11)
  S2p <- 1 - m2p            # post-bite survival in hut with LLIN protected people - Expression (14)
  
  
  ### Transition probability from HS (host-seeking) to F (fed) state and from F to HS
  Pd <- Du*Eu + Dp*Ep 			              # Proba that an HS vector will be diverted (to HS next night if it survive) - Expression (15)
  Sd <- S 	  									          # Proba that diverted vector survives to state HS next night - Expression (16)
  Pf_u <- Eu * fi_u								        # Proba that an HS vector will bite succesfully (same night) on an unprotected human - Expression (17)
  Pf_p <- Ep * fi_p			                  # Proba that an HS vector will bite succesfully (same night) on a LLIN protected human - Expression (18)
  Pf <- Pf_u + Pf_p										    # Proba that an HS vector will bite succesfully (same night) on human - Expression (19)
  F_u <- Pf_u / (Pf_u + Pf_p)							# Proportion fed on unprotected humans - Expression (20)
  Sf <- (S2u * F_u + (1-F_u)*S2p) * S^g		# Proba survives from successful bite to HS at next G cycle - Expression (21)
  
  
  ### model average lifetime infectious bites
  
  PfA <- Pf/(1-Pd*Sd)									  # Average proba that a HS vector will survive to take a feed - Expression (22)
  BA <- PfA/(1-Sf*PfA)							  	# Average number of bites which a HS vector will survive to give - Expression (23)
  c <- k * Ih          		              # probability that a vector become infected while taking a bloodmeal -  Expression (24)
  Pl <- PfA*c/(1-PfA*Sf*(1-c))					# Probability that a vector will acquire Pf during its lifetime - Expression (25)
  
  # Expression (35)

  move <- function(i, g, n){				# this fonction return 1 if the duration of the combination equals the duration to become infectious (n)
    if ((ceiling((n-g-i)/g)*g+i)==(n-g)){
      return (1)
    } else {
      return (0)
    }
  }
  
  Term1 <- NULL
  Term2 <- NULL
  Term3 <- NULL

  for (i in 0:(n-g)){
    Ng <- ceiling((n-g-i)/g)      # number of gonotrophic cycles required to complete sporogony (giving i diversion events)
    Ng2 <- ceiling((n-2*g-i)/g)   # number of gonotrophic cycles required to complete sporogony (giving i diversion events and the combination finshes by a g)
    Term1[i+1] <- ((Pf*Sf)^Ng) * ((Pd*Sd)^i)    # proba to survive the combination of i and g
    Term2[i+1] <- i+Ng2+move(i,g,n)     	      # n in the binomial coefficient (total number of element)
    Term3[i+1] <- factorial(Term2[i+1]) / (factorial(Term2[i+1]-i)*factorial(i))	# Binomial coefficient : nb of possible order of each combination (with 0 to n-g diversions events)
  }
  
  Term4 <- Term1*Term3	       # Term3 gives the number of possible orders for each combinaison of (i) diversions and ((n-g-i)/g) feeds
  Sl <- Sf * sum(Term4)				 # Proba that a newly infected vector will survive to HS state as an infectious vector - Expression (26)
  
  VLAIB = Pl*Sl*BA						 # Vector average lifetime infectious bites (= individual Vector capacity) - Expression (27)
  
  ### number of oviposition events - can be used as a fitness indicators in a genetic/evolution model - Expression (30)
  OvA <- (1/(1-(PfA*Sf)))-1    # Average number of oviposition which a HS vector will survive to give (geometric series of first term 1 and reason Sf*Pfa, decreased by 1)
  
  results <- c(Eu, Ep, Pd, Sd, Pf_u, Pf_p, Pf, F_u, Sf, Sl, PfA, BA, OvA, Pl, VLAIB)
  names(results) <- c("Eu", "Ep", "Pd", "Sd", "Pf_u", "Pf_p", "Pf", "F_u", "Sf", "Sl", "PfA", "BA","OvA", "Pl", "VLAIB")
  return(results)
}


# function to calculate RTP (Relative transmission potential)
#' `fRTP` return the relative transmission potetntial of a scenario against another taken as baseline. It is the ratio of VLAIB
#' between to scenarii. It can use various model that are given as `FUN` argument
#'
#' @param nsim number of simulation (used when FUN is model that used simulation)
#' @param m m1p: pre-bite feeding related mortality probability when faced to un LLIN protected host (default = 0.72)
#' @param m2 m1p at baseline
#' @param p Pllin: preference for LLIN protected humans (against unprotected humans) as recorded in a dual choice olfactometer (default = 0.5, i.e inert LLIN)
#' @param p2 Pllin at baseline 
#' @param D Diversion probability when entering a hut with LLIN (default = 0.3)
#' @param D2 D2 at baseline
#' @param Uh LLIN use rate in the population (default = 0.6)
#' @param Uh2 Uh at baseline (default = 0.6)
#' @param pi pi: proportion of exposure to bite that occurs during which LLIN is in use (default = 0.9)
#' @param pi2 pi at baseline (default = 0.9)
#' @param FUN the function used to calculate VLAIB in both scenarii, the average number of infectious bites that a vector gives according to a behavior and mortality model
#' @param fu fi1_u: successful feeding probability of pre-bite survivors (default = 0.55)
#'
#' @return the value of Relative transmission potential (RTP)
#' @export
#'
#' @examples
#' fRTP(Uh=0.8, Uh2=0.5)
#' 
fRTP <- function(nsim = 1000,m = 0.72, m2= 0.72, p = 0.5, p2 = 0.5, D = 0.3, D2= 0.3, Uh = 0.6, Uh2 = 0.6, pi = 0.9, pi2 = 0.9, FUN = VLAIB, fu = 0.55){
  RTP <- FUN(nsim, m1p = m , Pllin = p , Dp = D , Uh = Uh , pi = pi)["VLAIB"] / 
         FUN(nsim, m1p = m2, Pllin = p2, Dp = D2, Uh = Uh2, pi = pi2)["VLAIB"]
  names(RTP) <- "RTP"
  return(RTP)
}






HSBM <- function(nsim, S = 0.9, g = 3, m1u = 0.016, m1p = 0.5, m2u = 0.005, m2p = 0.21, fi1_u = 0.55, 
                 RR_fi1 = 0.56, NH = 1000, Ch = 0.6, pi = 0.9, Pllin = 0.5){
  
  ## proba of choice between unprotected and protected human
  Np <- Ch*pi*Nh                      # average number of people protected by an LLIN - Expression (1)
  Nu <- Nh - Np                        # average number of unprotected people - Expression (2)
  Cpp <- dhyper(2, round(Np), round(Nu), 2, log = FALSE) # - Expression (4)
  Cpu <- dhyper(1, round(Np), round(Nu), 2, log = FALSE) # - Expression (5)
  Cuu <- dhyper(0, round(Np), round(Nu), 2, log = FALSE) # - Expression (6)
  Eu <- Cuu + Cpu * (1 - Pllin)				  # proba to enter a house with an unprotected human - Expression (7)
  Ep <- Cpp + Cpu * Pllin			          # proba to enter a house with a protected human - Expression (8)
  
  
  ## survival and diversion probabilities
  S1u  <- 1 - m1u             # pre-bite survival in hut with unprotected people - Expression (10)
  S1p  <- 1 - m1p             # pre-bite survival in hut with LLIN protected people  - Expression (16)
  
  S2u <- 1 - m2u             # post-bite survival in hut with unprotected people - Expression (14)
  S2p <- 1  - m2p            # post-bite survival in hut with LLIN protected people - Expression (20)
  
  fi_u  <- S1u * fi1_u          # Successful feeding probability when entering a hut with an unprotected human - Expression (12)
  fi1_p <- fi1_u * RR_fi1       # - Expression (22)
  fi_p  <- S1p * fi1_p          # - Expression (18)
  
  D1u <- 1-fi1_u               # - Expression (11)
  D1p <- 1-fi1_u*RR_fi1         # - Expression (17)
  
  Du  <- S1u * D1u                  # probability of diversion (HS next day) when entering a hut with an unprotected human  - Expression (13)
  Dp  <- S1p * D1p                  # probability of postpone (HS next day) when entering a hut with a protected human  - Expression (19)
  
  
  ### Transition probability from HS (host-seeking) to F (fed) state and from F to HS
  Pd <- Du*Eu + Dp*Ep 			              # Proba that an HS vector will be diverted (to HS next night if it survive) - Expression (24)
  Sd <- S 	  									          # Proba that diverted vector survives to state HS next night - Expression (25)
  Pf_u <- Eu * fi_u								        # Proba that an HS vector will bite succesfully (same night) on an unprotected human - Expression (26)
  Pf_p <- Ep * fi_p			                  # Proba that an HS vector will bite succesfully (same night) on a LLIN protected human - Expression (27)
  Pf <- Pf_u + Pf_p										    # Proba that an HS vector will bite succesfully (same night) on human - Expression (28)
  F_u <- Pf_u / (Pf_u + Pf_p)							# Proportion fed on unprotected humans - Expression (29)
  Sf <- (S2u * F_u + (1-F_u)*S2p) * S^g		# Proba survives from successful bite to HS at next G cycle - Expression (30)
  
  return(c(Pd,Sd,Pf,Sf))
}