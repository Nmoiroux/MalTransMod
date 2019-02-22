#' `fitness_f_kdr` returns relative fitnesses based on the average number of oviposition events that females anopheles 
#' of the three genotypes for the kdr mutation are expected to performed during its lifetime
#'
#' @param Pref_kdr a vector of size three of which elements are the preference for protected LLIN users of the three genotypes (in the following order : RR, RS, SS)
#' @param m1p_kdr a vector of size three of which elements are the pre-bite mortality of the three genotypes (in the following order : RR, RS, SS) when entering a hut with a LLIN protected individual
#' @param m2p_kdr a vector of size three of which elements are the pre-bite mortality of the three genotypes (in the following order : RR, RS, SS) when entering a hut with an unprotected individual
#' @param FUN the function used to calculate the number of oviposition events (default = `VLAIB``)
#'
#' @return a vector of the relative fitnesses (calculated from the number of oviposition events) of the three genotypes (in the following order: RR RS SS)
#' @export
#'
#' @examples
#' Pllin_kdr <- c(0.6,0.5,0.5)        # from porciani et al 2016
#' m1p_kdr   <- c(0.05,0.5,0.95)      # from Diop et al 2015
#' m2p_kdr   <- c(0.005,0.005,0.005)  # default values
#' fitness_f_kdr(Pllin_kdr, m1p_kdr, m2p_kdr, VLAIB)

fitness_f_kdr <- function(Pllin_kdr, m1p_kdr, m2p_kdr, FUN=VLAIB){
 
  OvA_RR <- FUN( 100, m1p = m1p_kdr[1], m2p = m2p_kdr[1], Pllin = Pllin_kdr[1], RR_fi1 = RR_fi1_kdr[1], fi1_u = fi1_u_kdr[1] )["OvA"]
  OvA_RS <- FUN( 100, m1p = m1p_kdr[2], m2p = m2p_kdr[2], Pllin = Pllin_kdr[2], RR_fi1 = RR_fi1_kdr[2], fi1_u = fi1_u_kdr[2] )["OvA"]
  OvA_SS <- FUN( 100, m1p = m1p_kdr[3], m2p = m2p_kdr[3], Pllin = Pllin_kdr[3], RR_fi1 = RR_fi1_kdr[3], fi1_u = fi1_u_kdr[3] )["OvA"]
  
  OvA_kdr <- c(OvA_RR, OvA_RS, OvA_SS)    # number of oviposition that each genotype is expected to do
  W_F <- OvA_kdr / max(OvA_kdr)						# relative fitness
  return(W_F)
}

