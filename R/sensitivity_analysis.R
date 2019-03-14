### uncertainty and sensitivity analysis
library(pse)
library(tidyverse)
source("R/Fun_VLAIB_fRTP.R")

#### define range (min and max) values for each parameters ----
### from litterature
S_rg   <- list(min=0.61, max=0.98) # Silver 2008, Chapter 13
g_rg   <- list(min=2, max=6)				# Afrane et al. 2005, JME; Carnevale & Robert 2009
k_rg   <- list(min=0.02, max=0.2)  # Churcher et al. 2015, Nat. Comm
n_rg   <- list(min=8, max=16)			# Hien 2016, Plos Path; Ohm et al. 2018, Par.& Vec.
Uh_rg  <- list(min=0.2, max=0.8)		# Malaria Atlas Project 
Ih_rg  <- list(min=0.1, max=0.9)   # Malaria Atlas Project 
pi_rg  <- list(min=0.45, max=1)    # Cooke et al. 2015, Malaria Journal

### user defined
Nh_rg  <- list(min=300, max=5000)
   
### from Moiroux et al. 2017
## load summarised data of Experimental hut trials (EHT) from Moiroux et al. 2017
Data_moiroux <- read.delim("data/Data_moiroux.txt")	

Data_moiroux %>%
	dplyr::mutate(m1 = Tot_D_unfd / Total_unfd) %>% 	       # calculate pre-bite mortality
	dplyr::mutate(m2 = Tot_D_fed / Total_bfed) %>%  	       # calculate post-bite mortality
	dplyr::mutate(D = Tot_L_unfd / total) %>%   	           # calculate Diversion rate
	group_by(ITN) %>% 																			 # group by type of tretament (ITN, CTN or control)
	summarise_at(c("m1","m2","D"),funs(min,max), na.rm=TRUE) -> stat_moiroux # calculate min and max values of m1, m2 and D

## find min and max values of pre-bite, post-bite mortality and diversion for untreated nets and LLINs
stat_UTN <- stat_moiroux %>% filter(ITN=="no") %>% select(-ITN) %>% as.data.frame()
stat_ITN <- stat_moiroux %>% filter(ITN=="ITN")%>% select(-ITN) %>% as.data.frame()

m1u_rg <- list(min=stat_UTN[1,1], max=stat_UTN[1,4])
m2u_rg <- list(min=stat_UTN[1,2], max=stat_UTN[1,5])
m1p_rg <- list(min=stat_ITN[1,1], max=stat_ITN[1,4])
m2p_rg <- list(min=stat_ITN[1,2], max=stat_ITN[1,5])
Du_rg  <- list(min=stat_UTN[1,3], max=stat_UTN[1,6])
Dp_rg  <- list(min=stat_ITN[1,3], max=stat_ITN[1,6])

## set preference values for the comparison of transmission
pref <- 0.6					# in the sensitivity analysis, we will compare an attractive LLIN to
pref_ref <- 0.3     # a deterrent one

#### prepare functions and data for uncertainity analysis----

### modify fRTP function to be fed with all parameters of VLAIB function ----
fRTP_sens <- function(nsim, S = 0.9, g = 3, Nh = 1000, Ih = 0.5, k = 0.1, n = 11, m1u = 0.05,
                                                                                  m1p = 0.72, 
                                                                                  m2u = 0.005,
                                                                                  m2p = 0.21,
																																									Du = 0.43,
                                                                                  Dp = 0.3,
                                                                                  Uh = 0.6,
                                                                                  pi = 0.9, 
                                                                                  #Pllin = 0.5,  
                                                                                  FUN = VLAIB){
													  RTP <- FUN(nsim,S = S, g = g, Nh = Nh, Ih = Ih, k = k, n = n, 
																																			            m1u=m1u,
																																			            m1p=m1p, 
																																			            m2u=m2u,
																																			            m2p=m2p,
																																			  					Du=Du,
																																			            Dp=Dp,
																																			            Uh=Uh,
																																			            pi=pi, 
																																			            Pllin=pref)["VLAIB"] / +
													         FUN(nsim,S = S, g = g, Nh = Nh, Ih = Ih, k = k, n = n, 
																																			        		m1u=m1u,
																																			        		m1p=m1p, 
																																			        		m2u=m2u,
																																			        		m2p=m2p,
																																			        		Du=Du,
																																			        		Dp=Dp,
																																			        		Uh=Uh,
																																			        		pi=pi, 
																																			        		Pllin=pref_ref)["VLAIB"]
  return(-(1-RTP)*100)
}

### following tutorial of pse package ----
## names of the parameters
factors <- c("S", "g", "Nh", "Ih", "k", "n", 
             "m1u" ,    # 0 - 0.04 (EHT Moiroux)
             "m1p",    # 0.03 - 1 (EHT Moiroux)          # TESTED
             "m2u",     # 0 - 0.01 (EHT Moiroux)
             "m2p",     # 0 - 0.5 (EHT Moiroux)
             "Du",       # 0.15 - 0.8 (EHT Moiroux)
             "Dp",     # 0 - 1                          # TESTED
             "Uh",     # 0 - 1                          # TESTED
             "pi"     # 0 - 1                          # TESTED
             #"Pllin"   # 0.1 - 0.7                      # TESTED)
             ) 


## the probability density functions for each parameter
# discrete uniform probablity function (for parameters g and n)
qdunif<-function(p, min, max) floor(qunif(p, min, max))

q <- c("qunif",  #S
       "qdunif", #g
       "qunif",  #Nh
       "qunif",  #Ih
       "qunif",  #k
       "qdunif", #n
       "qunif",  #m1u
       "qunif", #m1p     # TESTED (fixed)
       "qunif",  #m2u
       "qunif",  #m2p
       "qunif", #Du
       "qunif", #RR_fi1 # TESTED (fixed)
       "qunif", #Uh     # TESTED (fixed)
       "qunif"#, #pi     # TESTED (fixed)
       #"qunif"  #Pllin
       )


## a list containing the lists with all the parameters to the density functions ----
q.arg <- list(S_rg, #S
              g_rg, #g
              Nh_rg, #Nh
              Ih_rg,  #Ih
              k_rg, #k
              n_rg,     #n
              m1u_rg,   #m1u
              m1p_rg,  #m1p # TESTED (fixed)
              m2u_rg,   #m2u
              m2p_rg, #m2p
							Du_rg, #Du
							Dp_rg,   #Dp # TESTED (fixed)
							Uh_rg,   #Uh # TESTED (fixed)
							pi_rg#,   #pi # TESTED (fixed)
              #list(min=0.2, max=0.4)  #Pllin
              )


# function modelRun encapsulates fRTP_sens function,in a manner to receive a data.frame containing ----
# all parameter combinations and returning the results in one array.
modelRun <- function (my.data) {
  return(mapply(fRTP_sens, 100, my.data[,1], my.data[,2], my.data[,3], my.data[,4], my.data[,5], my.data[,6],
                my.data[,7], my.data[,8], my.data[,9], my.data[,10]
                , my.data[,11], my.data[,12], my.data[,13], my.data[,14]#, my.data[,15]
                ))
}

# Generates the Latin Hypercube sampling for uncertainty and sensitivity analyses.----
myLHS <- pse::LHS(modelRun, factors, 500, q, q.arg, nboot=50)

# accessing the data and result data frames from an LHS object
pse::get.data(myLHS)

# plots the empirical cumulative density function----
pse::plotecdf(myLHS)

# produces a series of scatterplots from data----
pse::plotscatter(myLHS)

# plots the partial rank correlation coefficient from an LHS object----
pse::plotprcc(myLHS)

# Estimates the partial inclination coefficient of a model response in relation with all model input variables----
pse::pic(myLHS, nboot=40)

# In order to decide whether our sample size was adequate or insufficient, we calculate
# the Symmetric Blest Measure of Agreement (SBMA) between the PRCC coefficients
# of two runs with different sample sizes.
newLHS <- pse::LHS(modelRun, factors, 750, q, q.arg)
(mySbma <- pse::sbma(myLHS, newLHS)) 
# 0.8986226 # It is reasonable to expect agreements around 0.7 to 0.9 in well-behaved models, but two cases require attention.



# End ----
### define possible range of parameters from field data
Data_EHT <- read.delim("data/Data_EHT.txt")
Data_meta <- read.delim("data/DB_meta_AnGambiae.txt")
Data_EHT$Total_unfd <- Data_EHT$Tot_L_unfd + Data_EHT$Tot_D_unfd
Data_meta$Total_unfd <- Data_meta$Tot_L_unfd + Data_meta$Tot_D_unfd

# create a new column categorizing nets in three categorie (no=control, CTN and ITN)
Data_meta$ITN <- as.factor(Data_meta$ttmt2)
levels(Data_meta$ITN)[1] <- "no"
levels(Data_meta$ITN)[grep("CTN",levels(Data_meta$ITN))] <- "CTN"
levels(Data_meta$ITN)[!(levels(Data_meta$ITN) %in% c("no","CTN"))] <- "ITN"

# calculate m1, m2 and D 
Data_meta %>% dplyr::mutate(m1 = Tot_D_unfd / Total_unfd) -> Data_meta 
Data_meta %>% dplyr::mutate(m2 = Tot_D_fed / Total_bfed) -> Data_meta
Data_meta %>% dplyr::mutate(D = Tot_L_unfd / total) -> Data_meta

summarytools::view(summarytools::dfSummary(Data_meta[,c(2,22,29:32)] %>% dplyr::filter(ITN == "no", total > 0)))	
summarytools::view(summarytools::dfSummary(Data_meta[,c(2,22,29:32)] %>% dplyr::filter(ITN == "ITN", total > 0)))	

summarytools::view(summarytools::dfSummary(Data_EHT[,c(1:3,11:13)] %>% dplyr::filter(ITN == "no")))
summarytools::view(summarytools::dfSummary(Data_EHT[,c(1:3,11:13)] %>% dplyr::filter(ITN == "LLIN")))

Data_EHT %>% dplyr::filter(ITN == "LLIN") -> data_ITN
fit <- fitdist(data_ITN$Tot_D_unfd,"binom", start=list(size=data_ITN$Total_unfd, prob=data_ITN$m1))
fit <- fitdist(data_ITN$Tot_D_unfd,"binom", fix.arg=list(size=data_ITN$Total_unfd), start=list(prob=data_ITN$m1))

###### mise au propre
Data_EHT <- read.delim("data/Data_EHT.txt")
Data_EHT$Total_unfd <- Data_EHT$Tot_L_unfd + Data_EHT$Tot_D_unfd
Data_EHT %>% dplyr::filter(ITN == "no") %>% summarise_at(c("m1","m2","D"),funs(min,max))
Data_EHT %>% dplyr::filter(ITN == "LLIN") %>% summarise_at(c("m1","m2","D"),funs(min,max))









boxplot(subset(Data_EHT, ITN=="no")$m1)
boxplot(subset(Data_EHT, ITN=="LLIN")$m1)
boxplot(subset(Data_EHT, ITN=="no")$m2)
boxplot(subset(Data_EHT, ITN=="LLIN")$m2)
boxplot(subset(Data_EHT, ITN=="no")$fi1)
boxplot(subset(Data_EHT, ITN=="LLIN")$fi1)
boxplot(subset(Data_EHT, ITN=="LLIN")$RR_fi1)

hist(subset(Data_EHT, ITN=="no")$m1)
hist(subset(Data_EHT, ITN=="LLIN")$m1)
hist(subset(Data_EHT, ITN=="no")$m2)
hist(subset(Data_EHT, ITN=="LLIN")$m2)
hist(subset(Data_EHT, ITN=="no")$fi1)
hist(subset(Data_EHT, ITN=="LLIN")$fi1)


plot_EHT_LLIN <- Data_EHT[,c(1:3,11:14)] %>% gather("variable","value", 4:7) %>% filter(ITN == "LLIN")
ggplot2::ggplot(plot_EHT_LLIN, ggplot2::aes(x=variable, y=value)) + 
	ggplot2::geom_boxplot(alpha=0.4) +
	ggplot2::stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red")

plot_EHT_ctrl <- Data_EHT[,c(1:3,11:13)] %>% gather("variable","value", 4:6) %>% filter(ITN == "no", variable != "fi1")
ggplot2::ggplot(plot_EHT_ctrl, ggplot2::aes(x=variable, y=value)) + 
	ggplot2::geom_boxplot(alpha=0.4) +
	ggplot2::stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red")

plot_EHT_ctrl <- Data_EHT[,c(1:3,11:13)] %>% gather("variable","value", 4:6) %>% filter(ITN == "no", variable == "fi1")
ggplot2::ggplot(plot_EHT_ctrl, ggplot2::aes(x=variable, y=value)) + 
	ggplot2::geom_boxplot(alpha=0.4) +
	ggplot2::stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red")


summarytools::view(summarytools::dfSummary(Data_EHT[,c(1:3,11:13)] %>% filter(ITN == "no")))

