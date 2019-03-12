library(ggpubr)
library(tidyverse)
insertSource("R/my_ggarrange.R", package = "ggpubr")

# -----
p<- c(0.6,0.7,0.8)		      # preference for LLIN protected human (against unprotected human) 
Cov <- seq(0,1,by = 0.01)	  # coverage LLIN
m <- seq(0,1,by = 0.01)     # pre-bite feeding-related mortality
d <- seq(0,1,by = 0.01)     # diversion probability
Behav <- seq(0,1,by = 0.01) # proportion of exposure to bite during which ITN is in use
ref_pref <- 0.5							#
FUN <- VLAIB

#### figure 1 (coverage)
# panel A
p<- c(0.6,0.7)		      # preference for LLIN protected human (against unprotected human) 
cov <- seq(0,1,by = 0.01)	  # coverage LLIN
ref_pref <- 0.5							# preference value of the reference LLIN (inert = 0.5)
FUN <- VLAIB                # function used to calculate vectorial capacity

RTP.fit <- expand.grid(p.=p,cov.=cov)   # create table of all combinations of p and cov
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(cov.,p.),~fRTP(p=.y, Uh=.x, Uh2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig1_A <- ggplot(RTP.fit, aes(x=cov.*100, y=redVC)) + 
	xlab("LLIN Coverage (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(2,3),name="LLIN attraction")+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())


# panel B
p<- c(0.5,0.6,0.7)		      # preference for LLIN protected human (against unprotected human) 
cov <- seq(0,1,by = 0.01)	  # coverage LLIN
ref_pref <- 0.3							# preference value of the reference LLIN (inert = 0.5)
FUN <- VLAIB                # function used to calculate vectorial capacity

RTP.fit <- expand.grid(p.=p,cov.=cov)   # create table of all combinations of p and cov
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(cov.,p.),~fRTP(p=.y, Uh=.x, Uh2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig1_B <- ggplot(RTP.fit, aes(x=cov.*100, y=redVC)) + 
	xlab("LLIN Coverage (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(1,2,3),name="LLIN attraction")+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())

fig1_B

# arrange panels of figure 1
figure1 <- ggarrange(fig1_A + rremove("legend"), 
										 fig1_B + rremove("legend"), 
										 common.legend = TRUE,
										 plot_legend = 2,
										 labels = c("A", "B"),
										 ncol = 2, nrow = 1)

annotate_figure(figure1,
								bottom = text_grob("LLIN coverage (%)"),
								left = text_grob("% reduction in vectorial capacity", rot=90)#,
								#fig.lab = "Figure 1", fig.lab.face = "bold"
)
