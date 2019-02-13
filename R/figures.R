library(ggpubr)
library(tidyverse)
library(export)


p<- c(0.6,0.7,0.8)		      # preference for LLIN protected human (against unprotected human) 
Cov <- seq(0,1,by = 0.01)	  # coverage LLIN
m <- seq(0,1,by = 0.01)     # pre-bite mortality
f <- seq(0,1,by = 0.01)     # blood-succes rate ratio in pre-bite survivors (vs untreated)
Behav <- seq(0,1,by = 0.01) # proportion of exposure to bite during which ITN is in use
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
  xlab("LLIN Coverage (%)") + #ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(x))) +
  scale_linetype_discrete(name="vector preference \n for LLINs (vs. untreated nets)")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())#+
	#theme(axis.title.x = element_blank())

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
  xlab("Survival in dwelling (%)") + #ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(x))) +
  #scale_linetype_discrete(name="vector preference \n for LLINs (vs. untreated nets)")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())#+
  #theme(axis.title.x = element_blank())
return(fig)
}

fig_m_1 <- fig_m(0.5,p,m)
fig_m_2 <- fig_m(0.3,p,m)

#### Graph RTP vs pref(x) and f_h_rr_live 

fig_f <- function(ref_pref,a,b){
RTP.fit <- expand.grid(x=a,y=b)
for (i in 1:length(RTP.fit$x)){
  #RTP.fit$z1[i] <- fRTP(p=RTP.fit$x[i],f=RTP.fit$y[i],f2=RTP.fit$y[i],fu=0.55, p2=ref_pref, FUN=VC_bincoef)
  RTP.fit$z2[i] <- fRTP(p=RTP.fit$x[i],f=RTP.fit$y[i],f2=RTP.fit$y[i],fu=0.55, p2=ref_pref, FUN=FUN)
}
RTP.fit$redVC <- 1-RTP.fit$z2
fig <- ggplot(RTP.fit, aes(x=(1-y)*100, y=-redVC*100)) + 
  xlab("Diversion probability (%)") + #ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(x))) +
  #scale_linetype_discrete(name="vector preference \n for LLINs (vs. untreated nets)")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())#+
	#theme(axis.title.x = element_blank())
return(fig)
}

fig_f_1 <- fig_f(0.5,p,f)
fig_f_2 <- fig_f(0.3,p,f)

### Graph RTP vs pref(x) and Pii(y)

fig_pii <- function(ref_pref,a,b){
RTP.fit <- expand.grid(x=a,y=b)
for (i in 1:length(RTP.fit$x)){
  #RTP.fit$z1[i] <- fRTP(p=RTP.fit$x[i],pii=RTP.fit$y[i],pii2=RTP.fit$y[i],p2=ref_pref,FUN=VC_bincoef)
  RTP.fit$z2[i] <- fRTP(p=RTP.fit$x[i],pi=RTP.fit$y[i],pi2=RTP.fit$y[i],p2=ref_pref,FUN=FUN)
}


RTP.fit$redVC <- 1-RTP.fit$z2
fig <- ggplot(RTP.fit, aes(x=(1-y)*100, y=-redVC*100)) + 
  xlab("LLIN avoidance (%)") + #ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(x))) +
  #scale_linetype_discrete(name="vector preference \n for LLINs (vs. untreated nets)")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())#+
#theme(axis.title.x = element_blank())
return(fig)
}

fig_pii_1 <- fig_pii(0.5,p,Behav)
fig_pii_2 <- fig_pii(0.3,p,Behav)

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
					common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 2, nrow = 4)

ggarrange(fig_cov_1, 
					fig_cov_2, 
					fig_m_1, 
					fig_m_2, 
					fig_f_1,
					fig_f_2,
					fig_pii_1,
					fig_pii_2,
					common.legend = TRUE,
					labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
					ncol = 2, nrow = 4)

library(gridExtra)
grid.arrange(fig_cov_1+ rremove("legend"), fig_cov_2, nrow=1)

ggarrange(fig_cov_1 + rremove("legend"), 
					fig_cov_2 + rremove("legend") , nrow=1,ncol = 2)
					
# export all graphs in one ppt
graph2ppt(fig_cov_1, file = "Rplot", append = FALSE)
graph2ppt(fig_cov_2, file = "Rplot", append = TRUE)
graph2ppt(fig_m_1, file = "Rplot", append = TRUE)
graph2ppt(fig_m_2, file = "Rplot", append = TRUE)
graph2ppt(fig_f_1, file = "Rplot", append = TRUE)
graph2ppt(fig_f_2, file = "Rplot", append = TRUE)
graph2ppt(fig_pii_1, file = "Rplot", append = TRUE)
graph2ppt(fig_pii_2, file = "Rplot", append = TRUE)


