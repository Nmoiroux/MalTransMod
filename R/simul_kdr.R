####### study impact on Kdr frequency (with various)


Pllin_kdr    <- c(0.3,0.3,0.3)   # preference for LLIN protected of genotype RR, RS and RS respectively
m1p_kdr   <- c(0.05,0.5,0.95)    # pre-bite mortality when faced to an LLIN of genotype RR, RS and SS respectively (Diop et al 2015, permethrin)
m2p_kdr  <- c(0.005,0.005,0.005) # post-bite mortality when faced to an LLIN of genotype RR, RS and SS respectively 

success_not <- c(1, 1, 1)		      # taux de succès de passage à travers la moustiquaire non traitée (= 1 si pas de moustiquaire) (RR, RS and SS)
success_net <- c(0.5, 0.75, 0.5)	# taux de succès de passage à travers la moustiquaire traitée (RR, RS and SS) (Diop et al 2015,  permethrin)
biting_not <- c(0.55,0.55,0.55)		# taux de succès du repas de sang (sans traitement) (RR, RS and SS) (Diop et al, unpublished,  permethrin)
biting_net <- c(0.8,0.4,0.2)		  # taux de succès du repas de sang (suite à un contact avec insecticide) (RR, RS and SS) (Diop et al, unpublished,  permethrin)

N.gen <- 100



df_pref <- as.data.frame(matrix(rep(c(0.7,0.5,0.3),3), 3,3))
df_pref[4,] <- c(0.6,0.5,0.5)
df_kdr <- merge(data.frame(1:N.gen), df_pref, by=NULL)

v_kdr <- NULL
for (i in (1:nrow(df_pref))){
  W_F <- fitness_f_kdr(as.numeric(df_pref[i,]),m1p_kdr,m2p_kdr, success_not, success_net, biting_not, biting_net) # relative fitness
  v_kdr <- c(v_kdr, predict_kdr(W_F))
}
df_kdr$fkdr <- v_kdr
df_kdr$scn <- as.factor(paste(df_kdr$V1, df_kdr$V2, df_kdr$V3))
levels(df_kdr$scn) <- c("RR: 0.7, RS: 0.7, SS: 0.7", "RR: 0.6, RS: 0.5, SS: 0.5", "RR: 0.5, RS: 0.5, SS: 0.5", "RR: 0.3, RS: 0.3, SS: 0.3")

fig <- ggplot(df_kdr, aes(x=X1.N.gen-1, y=fkdr)) + 
    xlab("N generations") + ylab("Kdr allelic frequency") +
    geom_line(aes(linetype=df_kdr$scn)) +
    scale_linetype_discrete(name="vector kdr genotype preference \n for LLINs (vs. untreated nets)")+
    xlim(0,25) + #ylim(min(RTP.fit$redVC),-1)+
    theme(aspect.ratio=1) 
    #theme(axis.title.y = element_blank())#+
  #theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
