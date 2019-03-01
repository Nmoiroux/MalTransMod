library(ggpubr)
library(tidyverse)
library(export)

# -----
p<- c(0.6,0.7,0.8)		      # preference for LLIN protected human (against unprotected human) 
Cov <- seq(0,1,by = 0.01)	  # coverage LLIN
m <- seq(0,1,by = 0.01)     # pre-bite feeding-related mortality
d <- seq(0,1,by = 0.01)     # diversion probability
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
  xlab("LLIN Coverage (%)") + #ylab("% reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(x))) +
  scale_linetype_discrete(name="LLIN attraction")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())

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
  xlab("feeding attempt survival (%)") + #ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(x))) +
  scale_linetype_discrete(name="LLIN attraction")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())
return(fig)
}

fig_m_1 <- fig_m(0.5,p,m)
fig_m_2 <- fig_m(0.3,p,m)

#### Graph RTP vs pref(x) and f_h_rr_live 

fig_f <- function(ref_pref,a,b){
RTP.fit <- expand.grid(x=a,y=b)
for (i in 1:length(RTP.fit$x)){
  #RTP.fit$z1[i] <- fRTP(p=RTP.fit$x[i],f=RTP.fit$y[i],f2=RTP.fit$y[i],fu=0.55, p2=ref_pref, FUN=VC_bincoef)
  RTP.fit$z2[i] <- fRTP(p=RTP.fit$x[i],D=RTP.fit$y[i],D2=RTP.fit$y[i],fu=0.55, p2=ref_pref, FUN=FUN)
}
RTP.fit$redVC <- 1-RTP.fit$z2
fig <- ggplot(RTP.fit, aes(x=(1-y)*100, y=-redVC*100)) + 
  xlab("Escaping (%)") + #ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(x))) +
  scale_linetype_discrete(name="LLIN attraction")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())
return(fig)
}

fig_f_1 <- fig_f(0.5,p,d)
fig_f_2 <- fig_f(0.3,p,d)

### Graph RTP vs pref(x) and Pii(y)

fig_pii <- function(ref_pref,a,b){
RTP.fit <- expand.grid(x=a,y=b)
for (i in 1:length(RTP.fit$x)){
  #RTP.fit$z1[i] <- fRTP(p=RTP.fit$x[i],pii=RTP.fit$y[i],pii2=RTP.fit$y[i],p2=ref_pref,FUN=VC_bincoef)
  RTP.fit$z2[i] <- fRTP(p=RTP.fit$x[i],pi=RTP.fit$y[i],pi2=RTP.fit$y[i],p2=ref_pref,FUN=FUN)
}


RTP.fit$redVC <- 1-RTP.fit$z2
fig <- ggplot(RTP.fit, aes(x=(1-y)*100, y=-redVC*100)) + 
  xlab("Spatial-temporal avoidance (%)") + #ylab("fold-reduction in vectorial capacity") +
  geom_line(aes(linetype=as.factor(x))) +
  scale_linetype_discrete(name="LLIN attraction")+
  xlim(0,100) + ylim(-100,0)+
  theme(aspect.ratio=1) +
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())
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

figure1 <- ggarrange(fig_cov_1 + rremove("legend"), 
					fig_cov_2 + rremove("legend"), 
					common.legend = TRUE,
					labels = c("A", "B"),
					ncol = 2, nrow = 1)

annotate_figure(figure1,
								bottom = text_grob("LLIN coverage (%)"),
								left = text_grob("% reduction in vectorial capacity", rot=90)#,
								#fig.lab = "Figure 1", fig.lab.face = "bold"
)

fig_cov_1$data %>% group_by(x) %>% summarise(max=max(redVC)) -> maxVC
fig_cov_1$data %>% filter(redVC %in% maxVC$max)

fig_cov_2$data %>% group_by(x) %>% summarise(max=max(redVC)) -> maxVC
fig_cov_2$data %>% filter(redVC %in% maxVC$max)


figure2 <- ggarrange(fig_m_1 + rremove("legend"), 
										 fig_m_2 + rremove("legend"), 
										 common.legend = TRUE,
										 labels = c("A", "B"),
										 ncol = 2, nrow = 1)

annotate_figure(figure2,
								bottom = text_grob("Survival per feeding attempt (%)"),
								left = text_grob("% reduction in vectorial capacity", rot=90)#,
								#fig.lab = "Figure 1", fig.lab.face = "bold"
)

fig_m_1$data %>% group_by(x) %>% summarise(max=max(redVC)) -> maxVC
fig_m_1$data %>% filter(redVC %in% maxVC$max)

fig_m_2$data %>% group_by(x) %>% summarise(max=max(redVC)) -> maxVC
fig_m_2$data %>% filter(redVC %in% maxVC$max)


figure3 <- ggarrange(fig_f_1 + rremove("legend"), 
										 fig_f_2 + rremove("legend"), 
										 common.legend = TRUE,
										 labels = c("A", "B"),
										 ncol = 2, nrow = 1)

annotate_figure(figure3,
								bottom = text_grob("Escaping (%)"),
								left = text_grob("% reduction in vectorial capacity", rot=90)#,
								#fig.lab = "Figure 1", fig.lab.face = "bold"
)

fig_f_1$data %>% group_by(x) %>% summarise(max=max(redVC)) -> maxVC
fig_f_1$data %>% filter(redVC %in% maxVC$max)

fig_f_2$data %>% group_by(x) %>% summarise(max=max(redVC)) -> maxVC
fig_f_2$data %>% filter(redVC %in% maxVC$max)


figure4 <- ggarrange(fig_pii_1 + rremove("legend"), 
										 fig_pii_2 + rremove("legend"), 
										 common.legend = TRUE,
										 labels = c("A", "B"),
										 ncol = 2, nrow = 1)

annotate_figure(figure4,
								bottom = text_grob("Spatial-temporal avoidance (%)"),
								left = text_grob("% reduction in vectorial capacity", rot=90)#,
								#fig.lab = "Figure 1", fig.lab.face = "bold"
)

fig_f_1$data %>% filter(y == 0.5)

fig_f_2$data %>% filter(y == 0.5)


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

#####


####### study impact on Kdr frequency (with various)----

#### (better) parameters (from litterature) for fitness calculus in females per genotypes----
m1p_kdr   <- c(0.05,0.5,0.95)     # pre-bite mortality when faced to an LLIN of genotype RR, RS and SS respectively (Diop et al 2015, permethrin)
m2p_kdr   <- c(0.005,0.005,0.005) # post-bite mortality when faced to an LLIN of genotype RR, RS and SS respectively 

success_not <- c(1, 1, 1)		      # taux de succès de passage à travers la moustiquaire non traitée (= 1 si pas de moustiquaire) (RR, RS and SS)
success_net <- c(0.5, 0.75, 0.5)	# taux de succès de passage à travers la moustiquaire traitée (RR, RS and SS) (Diop et al 2015,  permethrin)
biting_not  <- c(0.55,0.55,0.55)	# taux de succès du repas de sang (sans traitement) (RR, RS and SS) (Diop et al, unpublished,  permethrin)
biting_net  <- c(0.8,0.4,0.2)		  # taux de succès du repas de sang (suite à un contact avec insecticide) (RR, RS and SS: 0.8,0.4,0.2) (Diop et al, unpublished,  permethrin)

fi1_p_kdr   <-  success_net*biting_net      # Successful feeding probability of alive mosquitoes in treatment (genotype RR, RS and SS respectively)
fi1_u_kdr   <-  success_not*biting_not      # Successful feeding probability of alive mosquitoes in control (genotype RR, RS and SS respectively)
RR_fi1_kdr  <-  fi1_p_kdr / fi1_u_kdr       # risk ratio of successful feeding for alive mosquitoes in a hut with LLIN (compared to hut without LLIN)

#### parameter for a most simple demonstration (only based on pre-bite mortality)
m1p_kdr     <- c(0.05,0.5,0.95)     # pre-bite mortality when faced to an LLIN of genotype RR, RS and SS respectively (Diop et al 2015, permethrin)
m2p_kdr     <- c(0.005,0.005,0.005) # post-bite mortality when faced to an LLIN of genotype RR, RS and SS respectively 
fi1_u_kdr   <- c(0.55, 0.55, 0.55)  # Successful feeding probability of alive mosquitoes in a hut without LLIN (genotype RR, RS and SS respectively)
RR_fi1_kdr  <- c(0.56, 0.56, 0.56)  # risk ratio of successful feeding for alive mosquitoes in a hut with LLIN (compared to hut without LLIN)
N.gen <- 100 												# nb of generations to simulate

# create a dataframe with different combinations of preference (Pllin) per genotypes
df_pref <- as.data.frame(matrix(rep(c(0.7,0.5,0.3),3), 3,3))  # same pref. value for all genotypes, various level of preference (attraction, neutral, repulsion)
df_pref[4,] <- c(0.6,0.5,0.5)																  # only RR genotype attracted (as in Poriciani et. al 2017)
# data frame that will receive the results
df_kdr <- merge(data.frame(1:N.gen), df_pref, by=NULL)				# 

# loop that calculate 
v_kdr <- NULL
for (i in (1:nrow(df_pref))){
	W_F <- fitness_f_kdr(as.numeric(df_pref[i,]),m1p_kdr,m2p_kdr) # relative fitness of females
	v_kdr <- c(v_kdr, predict_kdr4(W_F=W_F,W_M=c(1,1,1), N.gen=N.gen)) # vector of R allele frequency
	
}
df_kdr$fkdr <- v_kdr
df_kdr$scn <- as.factor(paste(df_kdr$V1, df_kdr$V2, df_kdr$V3))
levels(df_kdr$scn) <- c("deterrent for all genotypes", "inert for all genotypes", "attractive for RR and inert for others", "attractive for all genotypes")


fig <- ggplot(df_kdr, aes(x=X1.N.gen-1, y=fkdr)) + 
	xlab("Generation") + ylab("Kdr allelic frequency") +
	geom_line(aes(linetype=fct_rev(df_kdr$scn))) +
	scale_linetype_discrete(name="LLIN remote effect")+
	xlim(0,25) + #ylim(min(RTP.fit$redVC),-1)+
	theme(aspect.ratio=1) 
#theme(axis.title.y = element_blank())#+
#theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
fig

fig$data[1:15,] 
fig$data %>% filter(X1.N.gen %in% c(1:15))



