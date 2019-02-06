library(ggpubr)
library(tidyverse)
library(popbio)


p<- c(0.6,0.7,0.8)		      # preference for LLIN protected human (against unprotected human) 
Cov <- seq(0,1,by = 0.01)	  # coverage LLIN
m <- seq(0,1,by = 0.01)     # pre-bite mortality
f <- seq(0,1,by = 0.01)     # blood-succes rate ratio in pre-bite survivors (vs untreated)
Behav <- seq(0,1,by = 0.01) # 
ref_pref <- 0.5							#
FUN <- VLAIB

### graph lines RTP vs pref(x) and Ch(y) 

fig_cov <- function(ref_pref,a,b){
RTP.fit <- expand.grid(x=a,y=b)
for (i in 1:length(RTP.fit$x)){
  #RTP.fit$z1[i] <- fRTP(p=RTP.fit$x[i],Ch=RTP.fit$y[i],Ch2=RTP.fit$y[i],p2=ref_pref,FUN=VC_bincoef)
  RTP.fit$z2[i] <- fRTP(100000, p=RTP.fit$x[i],Uh=RTP.fit$y[i],Uh2=RTP.fit$y[i],p2=ref_pref,FUN=FUN)
}

RTP.fit$redVC <- 1-RTP.fit$z2
fig <- ggplot(RTP.fit, aes(x=y*100, y=-redVC*100)) + 
  xlab("LLIN Coverage (%)") + ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(x))) +
  scale_linetype_discrete(name="vector preference \n for LLINs (vs. untreated nets)")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())#+
  #theme(axis.title.x = element_text(size = rel(1.2), angle = 00))

return(fig)
}
fig_cov_1 <- fig_cov(0.5,p,Cov)
fig_cov_2 <- fig_cov(0.3,p,Cov)


### graph lines RTP vs pref(x) and mortality(y)

fig_m <- function(ref_pref,a,b){
RTP.fit <- expand.grid(x=a,y=b)
for (i in 1:length(RTP.fit$x)){
  #RTP.fit$z1[i] <- fRTP(p=RTP.fit$x[i],m=RTP.fit$y[i],m2=RTP.fit$y[i],p2=ref_pref,FUN=VC_bincoef)
  RTP.fit$z2[i] <- fRTP(p=RTP.fit$x[i],m=RTP.fit$y[i],m2=RTP.fit$y[i],p2=ref_pref,FUN=FUN)
}

RTP.fit$redVC <- 1-RTP.fit$z2
fig <- ggplot(RTP.fit, aes(x=(1-y)*100, y=-redVC*100)) + 
  xlab("Survival in dwelling (%)") + ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(1-x))) +
  scale_linetype_discrete(name="vector preference \n for LLINs (vs. untreated nets)")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())#+
  #theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
return(fig)
}

fig_m_1 <- fig_m(0.5,p,m)
fig_m_2 <- fig_m(0.7,p,m)

#### Graph RTP vs pref(x) and f_h_rr_live 

fig_f <- function(ref_pref,a,b){
RTP.fit <- expand.grid(x=a,y=b)
for (i in 1:length(RTP.fit$x)){
  #RTP.fit$z1[i] <- fRTP(p=RTP.fit$x[i],f=RTP.fit$y[i],f2=RTP.fit$y[i],fu=0.55, p2=ref_pref, FUN=VC_bincoef)
  RTP.fit$z2[i] <- fRTP(p=RTP.fit$x[i],f=RTP.fit$y[i],f2=RTP.fit$y[i],fu=0.55, p2=ref_pref, FUN=FUN)
}
RTP.fit$redVC <- 1-RTP.fit$z2
fig <- ggplot(RTP.fit, aes(x=(1-y)*100, y=-redVC*100)) + 
  xlab("LLIN additional avoidance") + ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(1-x))) +
  scale_linetype_discrete(name="vector preference \n for LLINs (vs. untreated nets)")+
  xlim(0,100) + ylim(min(RTP.fit$redVC),-1)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())#+
  #theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
return(fig)
}

fig_f_1 <- fig_f(0.5,p,f)
fig_f_2 <- fig_f(0.7,p,f)

### Graph RTP vs pref(x) and Pii(y)

fig_pii <- function(ref_pref,a,b){
RTP.fit <- expand.grid(x=a,y=b)
for (i in 1:length(RTP.fit$x)){
  #RTP.fit$z1[i] <- fRTP(p=RTP.fit$x[i],pii=RTP.fit$y[i],pii2=RTP.fit$y[i],p2=ref_pref,FUN=VC_bincoef)
  RTP.fit$z2[i] <- fRTP(p=RTP.fit$x[i],pii=RTP.fit$y[i],pii2=RTP.fit$y[i],p2=ref_pref,FUN=FUN)
}


RTP.fit$redVC <- -1/RTP.fit$z2
fig <- ggplot(RTP.fit, aes(x=(1-y)*100, y=redVC)) + 
  xlab("Exposure to bites when/where LLINs are not in use (%)") + ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(1-x))) +
  scale_linetype_discrete(name="vector preference \n for LLINs (vs. untreated nets)")+
  xlim(0,100) + ylim(min(RTP.fit$redVC),-1)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())#+
  #theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
return(fig)
}

fig_pii_1 <- fig_pii(0.5,p,Behav)
fig_pii_2 <- fig_pii(0.7,p,Behav)

fig_cov_1
fig_cov_2
fig_m_1
fig_m_2
fig_f_1
fig_f_2
fig_pii_1
fig_pii_2

ggarrange(fig_cov_1 + rremove("legend"), 
          fig_cov_2 + rremove("legend"), 
          fig_m_1+ rremove("legend"), 
          fig_m_2+ rremove("legend"), 
          fig_f_1+ rremove("legend"),
          fig_f_2+ rremove("legend"),
          fig_pii_1+ rremove("legend"),
          fig_pii_2+ rremove("legend"),
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 2, nrow = 4)

########################## graph vector life
S = 0.9
g = 3
m_pre_h_u = 0.016
m_pre_h_p = 0.5
m_post_h_u = 0.005
m_post_h_p = 0.21
p_m_pre = 0.9
f_h_u_live = 0.55
f_h_p_live = 0.31
Pop = 1000
prev = 0.5
k = 0.1
n = 11
HBI = 1
Ch = 0
pii = 0.9
pref_h_u = 0.5
S_pre_h_u  <- 1 - m_pre_h_u               # pre-bite survival in control
S_pre_h_p  <- 1 - m_pre_h_p               # pre-bite survival in treatment
S_h_u <- 1 - m_post_h_u                                   # proba post-bite survival (24h) after feeding in a hut with an unprotected host / baseline post-bite mortality (always close to one)
S_h_p <- 1  - m_post_h_p                                  # post-bite survival in treatment (= fed mosquitoes survival)
S_h_rr<- S_h_p / S_h_u                                    # LLIN post-bite survival relative risk (S_h_p = S_h_u * S_h_rr) (average EHT Moiroux 2017)
f_h_u  <- S_pre_h_u * f_h_u_live    # Successful feeding probability when entering a hut with an unprotected human (proportion of Fed in the control hut)(very variable: 0.2 to 0.8, mean 0.5)
f_h_p  <- S_pre_h_p * f_h_p_live
f_h_rr <- f_h_p / f_h_u				      # LLIN feeding relative risk (f_h_p = f_h_u * f_h_rr) (= 1-BFI)
pp_h_u  <- S_pre_h_u * (1-f_h_u_live) # probability of postpone (HS next day) when entering a hut with an unprotected human (proportion unfed in control hut)(very variable: 0.2 to 0.8)
pp_h_p  <- S_pre_h_p * (1-f_h_p_live) # probability of postpone (HS next day) when entering a hut with a protected human 
pp_h_rr <- pp_h_p / pp_h_u		        # LLIN postpone relative risk (pp_h_p = pp_h_u * pp_h_rr -> proportion unfed alive in treated hut)
c <- k * prev * HBI  		# probability that a vector become infectious (exposed and then infectious) while taking a bloodmeal
p_ent_h_u <- dhyper(2, round((1-Ch*pii)*Pop), round(Ch*pii*Pop), 2, log = FALSE) + 			# proba to enter a house with an unprotected human when looking for a human host
  dhyper(1, round((1-Ch*pii)*Pop), round(Ch*pii*Pop), 2, log = FALSE) * pref_h_u	
p_ent_h_p <- dhyper(0, round((1-Ch*pii)*Pop), round(Ch*pii*Pop), 2, log = FALSE) + 			# proba to enter a house with a protected human when looking for a human host
  dhyper(1, round((1-Ch*pii)*Pop), round(Ch*pii*Pop), 2, log = FALSE) * (1 - pref_h_u)	
Pd <- pp_h_u*(p_ent_h_u + p_ent_h_p * pp_h_rr)			# Proba that an HS vector will be postpone to HS next night.
Sd <- S #	  											# Proba that posponed vector survives to state HS next night
Pf_u <- p_ent_h_u * f_h_u								# Proba that an HS vector will bite succesfully (same night) on an unprotected human
Pf_p <- p_ent_h_p * f_h_u * f_h_rr						# Proba that an HS vector will bite succesfully (same night) on a LLIN protected human
Pf <- Pf_u + Pf_p										# Proba that an HS vector will bite succesfully (same night) on human.
F_u <- Pf_u / (Pf_u + Pf_p)								# Proportion fed on unprotected humans
Sf <- S_h_u * (F_u + S_h_rr * (1-F_u))* S^g				# Proba survives from successful bite to HS at next G cycle.
Sfm <- Sf^(1/g)                     # daily survival rate betwwen successful bite and stage HS at next G cycle
t_max=50
c <- 0.05
HStoHS <- Pd*Sd
HStoF<- Pf 
HStoD <- 1-(HStoHS+HStoF)
FtoD <- 1-Sfm^(g/2)

simul_nbites_rec <- function(t_max.=t_max,HStoHS.=HStoHS,HStoF.=HStoF,HStoD.=HStoD,FtoD.=FtoD,g.=g, c.=c, n.=n){
  
  vFto <- rbinom(t_max.,1,FtoD.)
  vHSto <- rmultinom(t_max., 1, c(HStoHS.,HStoF.,HStoD.))
  v<-NULL
  v[1] <- match(1,vHSto[,1])-1
  for (i in 2:t_max.){
    if (v[i-1] == 2){                # if vector is in D state (2)
      break
      
    } else if (v[i-1] == 0){         # if vector is in HS state (0)
      v[i] <- match(1,vHSto[,i])-1#HSto()#                   # binomial probability of F or D
      
    } else if (v[i-1] == 1){         # if vector is in F state (1)
      if (i>(g.-1)) {                            # if vector has the age to complete sporogony (as vector take one time step (day) to become F, only g-1 days as F is sufficient to complete gonotrophy)
        nF <- sum(v[(i-(g.-1)):(i-1)]==1) # count number of days as F in the past g days
        if (nF==(g.-1)) {                      # if gonotrophic cycle is completed
          v[i] <- vFto[i]*2#Fto(FtoD)*2#                  # binomial probability of D or HS
        } else {                                # if gonotrophic cycle is not completed
          v[i] <- vFto[i]+1#Fto(FtoD)+1#                  # binomial probability of D or F
        }
      } else {
        v[i] <- vFto[i]+1#Fto(FtoD)+1#                # binomial probability of D or F
      }
    }
  }
  
  p<-NULL # vector of infectious state of the mosquito (0:susceptible, 1:exposed or 2:infectious)
  vPar <- rbinom(t_max.,1,c.)
  if (v[1]==1){    # first night
    p[1] <- vPar[i]
    #Agents[1,3] <- rbinom(1,1,c)#
  } else {
    p[1] <- 0  
  }
  
  for (i in 2:t_max.){
    if (v[i-1] == 2){                                                      # if vector is in D state (2)
      break
      
    } else if (p[i-1] == 0 && v[i-1] == 0 && v[i] == 1){ # non-exposed vector become F
      #Agents[i,3] <- rbinom(1,1,c)           # binomial probability to become exposed
      p[i] <- vPar[i]
    } else if (p[i-1] == 1){                                           # vector was exposed
      if(i>(n.-1)){                                                                  # if vector is sufisciently aged to complete spororgony
        nS <- sum(p[(i-(n.-1)):(i-1)]==1)  # count number of days as Exposed in the past n days
        if (nS==(n.-1)) {                                              # if spororgony is completed
          p[i] <- 2
        } else {
          p[i] <- 1
        }
      } else {
        p[i] <- 1
      }
      
    } else if (p[i-1] == 2){               # vector was infectious
      p[i] <- 2
      
    } else {
      p[i] <- 0
    }  
  }
  b<-NULL
  b[1] <- 0
  for (i in 2:t_max.){
    if (v[i-1] == 2){                # if vector is in D state (2)
      break
      
    }else if (v[i-1]==0 && v[i]==1 && p[i]==2){
      b[i]<-1
    } else{
      b[i]<-0
    }
  }
  res <- data.frame(v=v,p=p,b=b)
  return(res)
}

s <-simul_nbites_rec()
s$time <- seq.int(nrow(s))
s$b[s$b==0]<-NA

fig <- ggplot(s[-nrow(s),], aes(x=time, y=v, colour=factor(p))) +
  geom_line(aes(group = 1))+
  geom_point(aes(x=time-1, y=b-1) )
fig





