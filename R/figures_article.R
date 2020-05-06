library(ggpubr)
library(tidyverse)
source("R/Fun_VLAIB_fRTP.R")
source("R/Fun_fitness_f_kdr.R")
source("R/Fun_predict_kdr.R")
insertSource("R/my_ggarrange.R", package = "ggpubr") # modified version of ggarrange that allow to take the legend of any arranged plots as common legend

# -----

FUN <- IC # function used to calculate vectorial capacity

#### figure 1 (coverage) ----
# common value for Fig1 to 4
legend_title <- "Vector preference for LLINs"
y_title <- "% reduction in vectorial capacity"
# common parameters
cov <- seq(0,1,by = 0.01)	  # coverage LLIN
FUN <- IC                # function used to calculate vectorial capacity
# panel A
p<- c(0.6,0.7)		      # preference for LLIN protected human (against unprotected human) 
ref_pref <- 0.5							# preference value of the reference LLIN (inert = 0.5)


RTP.fit <- expand.grid(p.=p,cov.=cov)   # create table of all combinations of p and cov
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(cov.,p.),~fRTP(p=.y, Uh=.x, Uh2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig1_A <- ggplot(RTP.fit, aes(x=cov.*100, y=redVC)) + 
	xlab("LLIN Coverage (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(2,3),name=legend_title)+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())

# panel B
p<- c(0.5,0.6,0.7)		      # preference for LLIN protected human (against unprotected human) 
ref_pref <- 0.36							# preference value of the reference LLIN (attractive = 0.3)

RTP.fit <- expand.grid(p.=p,cov.=cov)   # create table of all combinations of p and cov
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(cov.,p.),~fRTP(p=.y, Uh=.x, Uh2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig1_B <- ggplot(RTP.fit, aes(x=cov.*100, y=redVC)) + 
	xlab("LLIN Coverage (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(1,2,3),name=legend_title)+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())


# arrange panels of figure 1
figure1 <- ggarrange(fig1_A + rremove("legend"), 
										 fig1_B + rremove("legend"), 
										 common.legend = TRUE,
										 plot_legend = 2,
										 labels = c("A", "B"),
										 ncol = 2, nrow = 1)

figure1 <- annotate_figure(figure1,
								bottom = text_grob("LLIN coverage (%)"),
								left = text_grob(y_title, rot=90))

#### figure 2 (physiological resistance) ----
# common parameters
m <- seq(0,1,by = 0.01)	  # pre-bite mortality

# panel A
p<- c(0.6,0.7)		      # preference for LLIN protected human (against unprotected human) 
ref_pref <- 0.5							# preference value of the reference LLIN (inert = 0.5)


RTP.fit <- expand.grid(p.=p,m.=m)   # create table of all combinations of p and m
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(m.,p.),~fRTP(p=.y, m=.x, m2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig2_A <- ggplot(RTP.fit, aes(x=(1-m.)*100, y=redVC)) + 
	xlab("Feeding attempt survival (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(2,3),name=legend_title)+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())

# panel B
p<- c(0.5,0.6,0.7)		      # preference for LLIN protected human (against unprotected human) 
ref_pref <- 0.36							# preference value of the reference LLIN (inert = 0.5)

RTP.fit <- expand.grid(p.=p,m.=m)   # create table of all combinations of p and m
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(m.,p.),~fRTP(p=.y, m=.x, m2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig2_B <- ggplot(RTP.fit, aes(x=(1-m.)*100, y=redVC)) + 
	xlab("feeding attempt survival (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(1,2,3),name=legend_title)+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())


# arrange panels of figure 1
figure2 <- ggarrange(fig2_A + rremove("legend"), 
										 fig2_B + rremove("legend"), 
										 common.legend = TRUE,
										 plot_legend = 2,
										 labels = c("A", "B"),
										 ncol = 2, nrow = 1)

figure2 <- annotate_figure(figure2,
								bottom = text_grob("Feeding attempt survival (%)"),
								left = text_grob(y_title, rot=90))

#### figure 3 (escaping) ----
d <- seq(0,1,by = 0.01)	  # diversion probability
# panel A
p<- c(0.6,0.7)		        # preference for LLIN protected human (against unprotected human) 
ref_pref <- 0.5							# preference value of the reference LLIN (inert = 0.5)

RTP.fit <- expand.grid(p.=p,d.=d)   # create table of all combinations of p and d
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(d.,p.),~fRTP(p=.y, D=.x, D2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig3_A <- ggplot(RTP.fit, aes(x=d.*100, y=redVC)) + 
	xlab("Escaping (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(2,3),name=legend_title)+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())

# panel B
p<- c(0.5,0.6,0.7)		      # preference for LLIN protected human (against unprotected human) 
ref_pref <- 0.36							# preference value of the reference LLIN (inert = 0.5)

RTP.fit <- expand.grid(p.=p,d.=d)   # create table of all combinations of p and d
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(d.,p.),~fRTP(p=.y, D=.x, D2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig3_B <- ggplot(RTP.fit, aes(x=d.*100, y=redVC)) + 
	xlab("Escaping (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(1,2,3),name=legend_title)+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())


# arrange panels of figure 1
figure3 <- ggarrange(fig3_A + rremove("legend"), 
										 fig3_B + rremove("legend"), 
										 common.legend = TRUE,
										 plot_legend = 2,
										 labels = c("A", "B"),
										 ncol = 2, nrow = 1)

figure3 <- annotate_figure(figure3,
								bottom = text_grob("Escaping (%)"),
								left = text_grob(y_title, rot=90))

#### figure 4 (spatial-temporal avoidance) ----
pi <- seq(0,1,by = 0.01)	  # proportion of exposure to bite during which LLIN is in use
# panel A
p<- c(0.6,0.7)		        # preference for LLIN protected human (against unprotected human) 
ref_pref <- 0.5							# preference value of the reference LLIN (inert = 0.5)

RTP.fit <- expand.grid(p.=p,pi.=pi)   # create table of all combinations of p and pi
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(pi.,p.),~fRTP(p=.y, pi=.x, pi2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig4_A <- ggplot(RTP.fit, aes(x=(1-pi.)*100, y=redVC)) + 
	xlab("Spatial-temporal avoidance (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(2,3),name=legend_title)+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())

# panel B
p<- c(0.5,0.6,0.7)		      # preference for LLIN protected human (against unprotected human) 
pi <- seq(0,1,by = 0.01)	  # 
ref_pref <- 0.36							# preference value of the reference LLIN (inert = 0.5)
FUN <- IC                # function used to calculate vectorial capacity

RTP.fit <- expand.grid(p.=p,pi.=pi)   # create table of all combinations of p and d
RTP.fit <- mutate(RTP.fit, z = pmap_dbl(list(pi.,p.),~fRTP(p=.y, pi=.x, pi2=.x, p2=ref_pref,FUN= FUN))) # calculate VC ratio
RTP.fit$redVC <- -(1-RTP.fit$z)*100        # reduction in VC

fig4_B <- ggplot(RTP.fit, aes(x=(1-pi.)*100, y=redVC)) + 
	xlab("Spatial-temporal avoidance (%)") + 
	geom_line(aes(linetype=as.factor(p.))) +
	scale_linetype_manual(values=c(1,2,3),name=legend_title)+
	xlim(0,100) + ylim(-100,0)+
	theme(aspect.ratio=1) +
	theme(axis.title.y = element_blank())+
	theme(axis.title.x = element_blank())


# arrange panels of figure 4
figure4 <- ggarrange(fig4_A + rremove("legend"), 
										 fig4_B + rremove("legend"), 
										 common.legend = TRUE,
										 plot_legend = 2,
										 labels = c("A", "B"),
										 ncol = 2, nrow = 1)

figure4 <- annotate_figure(figure4,
								bottom = text_grob("Spatial-temporal avoidance (%)"),
								left = text_grob(y_title, rot=90))

#### figure 5 (kdr dynamics) ----
#### parameters (only based on pre-bite mortality)
m1p_kdr     <- c(0.05,0.5,0.95)     # pre-bite mortality when faced to an LLIN of genotype RR, RS and SS respectively (Diop et al 2015, permethrin)
m2p_kdr     <- c(0.005,0.005,0.005) # post-bite mortality when faced to an LLIN of genotype RR, RS and SS respectively 
N.gen <- 100 												# nb of generations to simulate

# create a dataframe with different combinations of preference (Pllin) per genotypes
df_pref <- as.data.frame(matrix(rep(c(0.6,0.5,0.36),3), 3,3))  # same pref. value for all genotypes, various level of preference (attraction, neutral, deterrence)
df_pref[4,] <- c(0.6,0.5,0.5)																  # only RR genotype attracted (as in Poriciani et. al 2017)

# loop that calculate relative fitness for each combination of preference and simulate kdr dynamic accordingly
v_kdr <- NULL                                                        # vector that will receive allelic frequency values
for (i in (1:nrow(df_pref))){
	W_F <- fitness_f_kdr(as.numeric(df_pref[i,]),m1p_kdr,m2p_kdr)      # relative fitness of females
	v_kdr <- c(v_kdr, predict_kdr4(W_F=W_F,W_M=c(1,1,1), N.gen=N.gen)) # vector of R allele frequency
}

# data frame that will receive the results
df_kdr <- merge(data.frame(1:N.gen), df_pref, by=NULL)				
df_kdr$fkdr <- v_kdr

# prepare data for the plot's legend
df_kdr$scn <- as.factor(paste(df_kdr$V1, df_kdr$V2, df_kdr$V3))
levels(df_kdr$scn) <- c("deterrent for all genotypes", "inert for all genotypes", "attractive for RR and inert for others", "attractive for all genotypes")

# plot 
figure5 <- ggplot(df_kdr, aes(x=X1.N.gen-1, y=fkdr)) + 
	xlab("Number of generations") + ylab("Kdr allelic frequency") +
	geom_line(aes(linetype=fct_rev(df_kdr$scn))) +
	scale_linetype_discrete(name="LLIN remote effect")+
	xlim(0,25) + 
	theme(aspect.ratio=1) 




#### extraction of values of interest from the simulations (fig 1 to 5) ----

## search coverage value for which reduction in transmission is the higher
fig1_A$data %>% group_by(p.) %>% summarise(min=min(redVC)) -> minVC
fig1_A$data %>% filter(redVC %in% minVC$min)

fig1_B$data %>% group_by(p.) %>% summarise(min=min(redVC)) -> minVC
fig1_B$data %>% filter(redVC %in% minVC$min)

## search reduction in transmission value for survival set to 0 (no resistance)
fig2_A$data %>% filter(m. == 1)
fig2_B$data %>% filter(m. == 1)
## search reduction in transmission value for survival set to 1 (max resistance)
fig2_A$data %>% filter(m. == 0)
fig2_B$data %>% filter(m. == 0)

## search reduction in transmission value for escaping set to 0 (no quantitative bahavioral resistance)
fig3_A$data %>% filter(d. == 0)
fig3_B$data %>% filter(d. == 0)
## search reduction in transmission value for escaping set to 1 (max quantitative bahavioral resistance)
fig3_A$data %>% filter(d. == 1)
fig3_B$data %>% filter(d. == 1)

## search reduction in transmission value for spatial-temporal avoidance set to 0 (no qualitative bahavioral resistance)
fig4_A$data %>% filter(pi. == 1)
fig4_B$data %>% filter(pi. == 1)
## search reduction in transmission value for spatial-temporal avoidance set to 0.5 (50 % of exposure to bite occurs when LLIN are not in use)
fig4_A$data %>% filter(pi. == 0.5)
fig4_B$data %>% filter(pi. == 0.5)

## search number of generations needed to reach allelic frequency of 0.8
figure5$data %>% group_by(scn) %>% summarise(val = which.min(abs(fkdr-0.8)))


#### Supplementary Figure 1 plot preference value recorded by EHT in nature from two meta-analysis studies ----
# calculate prefernce value from field data
# from EHT in Moiroux et al. 2017, Plos One
Data_moiroux <- read.delim("data/Data_moiroux.txt")
Data_moiroux %>% 
	dplyr::filter(ITN == "no") %>% 
	select(Eval, total) %>% 
	right_join(Data_moiroux,by="Eval") %>% 
	#dplyr::filter(ITN == "ITN") %>%
	dplyr::filter(ITN != "no") %>%
	select(ttmt,total.x,total.y) %>%
	mutate(pLLIN = total.y/(total.y+total.x)) %>%
	mutate(study = "moiroux") %>%
	rename(Net = ttmt, TotalUTN = total.x, TotalITN = total.y) -> Data_pref_moiroux

# from EHT in Strode et al. 2014, Plos Med
Data_strode <- read.delim("data/Data_strode.txt")
Data_strode %>% 
	#filter(Intervention=="LLIN") %>%
	select(Net, TotalITN, TotalUTN) %>%
	mutate(study = "strode") %>%
	mutate(pLLIN = TotalITN / (TotalITN+TotalUTN)) -> Data_pref_strode

# bind the two dataframe together											 
Data_pref <- rbind(Data_pref_moiroux, Data_pref_strode)

# plot
supp_fig1 <- qplot( x=study , y=pLLIN , data=Data_pref , geom=c("boxplot","jitter") ,  ylab="Preference for treated net (pLLIN)")

# determine reference value of preference for deterrent LLIN (mean of pref among ITN that were significantly deterrent in Moiroux et al. 2017)
Data_deterence <- read.delim("data/Data_Figures2_Moiroux.txt")
Data_deterence %>% 
	filter(IC_high < 1) %>%	       # select only significantly deterent ITNs
	mutate(pLLIN=RR/(RR+1)) %>%    # convert Rate-ratio to preference values
	summarise(mean= mean(pLLIN)) %>%
	round(2)-> pLLIN_ref_det	 # compute the mean preference value

#### Supplementary Figure 2 plot Diversion, pre-bite and post-bite mortality value from Moiroux et al. 2017 ----
#### justify baseline value (m1, m2, D) from moiroux 2017 :
Data_moiroux <- read.delim("data/Data_moiroux.txt")	

Data_moiroux %>%
	dplyr::mutate(m1 = Tot_D_unfd / Total_unfd) %>% 	       # calculate pre-bite mortality
	dplyr::mutate(m2 = Tot_D_fed / Total_bfed) %>%  	       # calculate post-bite mortality
	dplyr::mutate(D = Tot_L_unfd / total) -> Data_moiroux2   # calculate Diversion rate

Data_moiroux2 %>%   	           
	group_by(ITN) %>% 																# group by type of tretament (ITN, CTN or control)
	summarise_at(c("m1","m2","D"),median, na.rm=TRUE) # calculate means of m1, m2 and D

# prepare data for ploting
Data_moiroux2 %>% 
	select(ttmt, ITN, m1, m2, D) %>%
	gather("m1","m2","D", key ="parameter", value = "value") -> plot_parameters

plot_parameters %>%
	filter(ITN == "no") -> plot_param_UTN

plot_parameters %>%
	filter(ITN == "ITN") -> plot_param_ITN

#plot
fig_param_A <- ggplot(plot_param_UTN, aes(x=parameter, y=value, fill=parameter)) +
	geom_boxplot(alpha=0.4) +
	geom_jitter()+
	stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
	theme(legend.position="none") +
	scale_fill_brewer(palette="Set3") +
	labs(y=" ")+
	ylim(0,1)+
	scale_x_discrete(labels=c(expression(paste('D'['U'])), expression(paste('µ'['1U'])), expression(paste('µ'['2U']))))+ 
	theme(axis.title.x = element_blank())

fig_param_B <- ggplot(plot_param_ITN, aes(x=parameter, y=value, fill=parameter)) +
	geom_boxplot(alpha=0.4)  +
	geom_jitter()+
	stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
	theme(legend.position="none") +
	scale_fill_brewer(palette="Set3") +
	labs(y=" ")+
	ylim(0,1)+
	scale_x_discrete(labels=c(expression(paste('D'['P'])), expression(paste('µ'['1P'])), expression(paste('µ'['2P']))))+
	theme(axis.title.x = element_blank())

figure_param <- ggarrange(fig_param_A, 
													fig_param_B, 
													labels = c("A", "B"),
													ncol = 2, nrow = 1)

annotate_figure(figure_param,
								bottom = text_grob("Parameters"),
								left = text_grob("value", rot=90)#,
								#fig.lab = "Figure 1", fig.lab.face = "bold"
)
