#' `fitness_f_kdr` returns relative fitness based on the average number of oviposition events that females anopheles 
#' of the three genotypes for the kdr mutation is expected to performed during its lifetime
#'
#' @param Pref_kdr a vector of size three of which elements are the preference for protected LLIN users of the three genotypes (in the following order : RR, RS, SS)
#' @param m1p_kdr a vector of size three of which elements are the pre-bite mortality of the three genotypes (in the following order : RR, RS, SS) when entering a hut with a LLIN protected individual
#' @param m2p_kdr a vector of size three of which elements are the pre-bite mortality of the three genotypes (in the following order : RR, RS, SS) when entering a hut with an unprotected individual
#' @param success_not a vector of size three of which elements are the success rate at penetrating through an holed untreated net (=1 if no net)
#' @param success_net a vector of size three of which elements are the success rate at penetrating through an holed LLIN
#' @param biting_not a vector of size three of which elements are the blood feeding success rate of pre-bite survivors in a hut with no (or untreated) net
#' @param biting_net a vector of size three of which elements are the blood feeding success rate of pre-bite survivors in a hut with an LLIN
#' @param FUN the function used to calculate the number of oviposition events (default = `VLAIB``)
#'
#' @return a vector of the relative fitnesses (calculated from the number of oviposition events) of the three genotypes (in the following order: RR RS SS)
#' @export
#'
#' @examples
#' Pllin_kdr <- c(0.6,0.5,0.5)        # from porciani et al 2016
#' m1p_kdr   <- c(0.05,0.5,0.95)      # from Diop et al 2015
#' m2p_kdr   <- c(0.005,0.005,0.005)  # default values
#' success_not <- c(1, 1, 1)		      # 
#' success_net <- c(0.5, 0.75, 0.5)	  # from Diop et al 2015
#' biting_not <- c(0.55,0.55,0.55)		# from Diop et al, unpublished
#' biting_net <- c(0.8,0.4,0.2)		    # from Diop et al, unpublished
#' fitness_f_kdr(Pllin_kdr, m1p_kdr, m2p_kdr, success_not, success_net, biting_not, biting_net, VLAIB)

fitness_f_kdr <- function(Pllin_kdr, m1p_kdr, m2p_kdr, success_not, success_net, biting_not, biting_net, FUN=VLAIB){
  
  fi1_p_kdr   <-  success_net*biting_net      # Successful feeding probability of alive mosquitoes in treatment (genotype RR, RS and SS respectively)
  fi1_u_kdr   <-  success_not*biting_not      # Successful feeding probability of alive mosquitoes in control (genotype RR, RS and SS respectively)
  RR_fi1_kdr  <-  fi1_p_kdr / fi1_u_kdr       # risk ratio of successful feeding for alive mosquitoes in a hut with LLIN (compared to hut without LLIN)
  
  OvA_RR <- FUN( 100, m1p = m1p_kdr[1], m2p = m2p_kdr[1], Pllin = Pllin_kdr[1], RR_fi1 = RR_fi1_kdr[1], fi1_u = fi1_u_kdr[1] )["OvA"]
  OvA_RS <- FUN( 100, m1p = m1p_kdr[2], m2p = m2p_kdr[2], Pllin = Pllin_kdr[2], RR_fi1 = RR_fi1_kdr[2], fi1_u = fi1_u_kdr[2] )["OvA"]
  OvA_SS <- FUN( 100, m1p = m1p_kdr[3], m2p = m2p_kdr[3], Pllin = Pllin_kdr[3], RR_fi1 = RR_fi1_kdr[3], fi1_u = fi1_u_kdr[3] )["OvA"]
  
  OvA_kdr <- c(OvA_RR, OvA_RS, OvA_SS)    # number of oviposition that each genotype is expected to do
  W_F <- OvA_kdr / max(OvA_kdr)						# relative fitness
  return(W_F)
}

