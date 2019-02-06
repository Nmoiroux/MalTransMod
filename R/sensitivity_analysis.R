### sensitivity analysis
library(pse)



### modify fRTP function to be fed with parameters of VC_bincoef function
fRTP_sens <- function(nsim,S = 0.9, g = 3, Pop = 1000, prev = 0.5, k = 0.1, n = 11, 
                                                                                    m_pre_h_u,
                                                                                    #m_pre_h_p, 
                                                                                    m_post_h_u,
                                                                                    m_post_h_p,
                                                                                    f_h_u_live,
                                                                                    #f_h_rr_live,
                                                                                    #Ch,
                                                                                    #pii, 
                                                                                    #pref_h_u,  
                                                                                    FUN){
  RTP<- FUN(nsim,S = S, g = g, Pop = Pop, prev = prev, k = k, n = n, 
            m_pre_h_u=m_pre_h_u,
            #m_pre_h_p=m_pre_h_p, 
            m_post_h_u= m_post_h_u,
            m_post_h_p=m_post_h_p,
            f_h_u_live=f_h_u_live,
            #f_h_rr_live=f_h_rr_live,
            #Ch=Ch,
            #pii=pii, 
            pref_h_u=0.3
            )["VLAIB",] / +
        FUN(nsim,S = S, g = g, Pop = Pop, prev = prev, k = k, n = n, 
            m_pre_h_u=m_pre_h_u,
            #m_pre_h_p=m_pre_h_p, 
            m_post_h_u= m_post_h_u,
            m_post_h_p=m_post_h_p,
            f_h_u_live=f_h_u_live,
            #f_h_rr_live=f_h_rr_live,
            #Ch=Ch,
            #pii=pii, 
            pref_h_u=0.5
            )["VLAIB",]
  return(RTP)
}

##### following tutorial of pse package
## names of the parameters
factors <- c("S", "g", "Pop", "prev", "k", "n", 
             "m_pre_h_u" , # 0 - 0.04 (EHT Moiroux)
             #"m_pre_h_p",  # 0.1 - 1 (EHT Moiroux)           # TESTED
             "m_post_h_u", # 0 - 0.01 (EHT Moiroux)
             "m_post_h_p", # 0 - 0.5 (EHT Moiroux)
             "f_h_u_live" # 0.2 - 0.8 (EHT Moiroux)
             #"f_h_rr_live", # 0.25 - 0.9 (EHT Moiroux)       # TESTED
             #"Ch", # 0 - 1                                  # TESTED
             #"pii", # 0 - 1                                 # TESTED
             #"pref_h_u" # 0.3 - 0.9                        # TESTED)
             ) 


## the probability density functions
# discrete uniform probablity function (for parameters g and n)
qdunif<-function(p, min, max) floor(qunif(p, min, max))

q <- c("qunif", #S
       "qdunif",#g
       "qunif", #pop
       "qunif", #prev
       "qunif", #k
       "qdunif", #n
       "qunif", #m_pre_h_u
       #"qunif", #m_pre_h_p # TESTED (fixed)
       "qunif", # m_post_h_u
       "qunif", # m_post_h_p
       "qunif"#, #f_h_u_live
       #"qunif", # f_h_rr_live # TESTED (fixed)
       #"qunif", # Ch # TESTED (fixed)
       #"qunif", # pii # TESTED (fixed)
       #"qunif"  #pref_h_u
       )
## a list containing the lists with all the parameters to the density functions
q.arg <- list(list(min=0.7, max=0.95), #S
              list(min=2, max=6),      #g
              list(min=300, max=5000), #pop
              list(min=0.1, max=0.8),  #prev
              list(min=0.02, max=0.2), #k
              list(min=8, max=16),     #n
              list(min=0, max=0.04),   #m_pre_h_u
              #list(min=0.1, max=0.9), #m_pre_h_p # TESTED (fixed)
              list(min=0, max=0.02),  # m_post_h_u
              list(min=0.08, max=0.5),# m_post_h_p
              list(min=0.2, max=0.8)#,  #f_h_u_live
              #list(min=0.25, max=0.9),  # f_h_rr_live # TESTED (fixed)
              #list(min=0.01, max=1),    # Ch # TESTED (fixed)
              #list(min=0.1, max=1),    # pii # TESTED (fixed)
              #list(min=0.2, max=0.4)   #pref_h_u
              )

# function that calculate de reduction in RTP
oneRun <- function(S, g, Pop, prev, k, n,
                   m_pre_h_u,
                   #m_pre_h_p, 
                   m_post_h_u,
                   m_post_h_p,
                   f_h_u_live#,
                   #f_h_rr_live,
                   #Ch,
                   #pii, 
                   #pref_h_u
                   ){
  x <- -1/ fRTP_sens(nsim,S, g, Pop, prev, k, n, 
                     m_pre_h_u,
                     #m_pre_h_p, 
                     m_post_h_u,
                     m_post_h_p,
                     f_h_u_live#,
                     #f_h_rr_live,
                     #Ch,
                     #pii, 
                     #pref_h_u
                     , FUN)
  return(x)
  }

# applicate oneRun on a dataframe
modelRun <- function (my.data) {
  return(mapply(oneRun, my.data[,1], my.data[,2], my.data[,3], my.data[,4], my.data[,5], my.data[,6],
                my.data[,7], my.data[,8], my.data[,9], my.data[,10]
                #, my.data[,11], my.data[,12], my.data[,13], my.data[,14], my.data[,15]
                ))
  }
myLHS <- LHS(modelRun, factors, 500, q, q.arg, nboot=50)
get.data(myLHS)
plotecdf(myLHS)
plotscatter(myLHS)
plotprcc(myLHS)
pic(myLHS, nboot=40)
newLHS <- LHS(modelRun, factors, 750, q, q.arg)
(mySbma <- sbma(myLHS, newLHS))

### define possible range of parameters

boxplot(subset(Data_EHT, ITN=="no")$m_pre_h_u.p)
boxplot(subset(Data_EHT, ITN=="LLIN")$m_pre_h_u.p)
boxplot(subset(Data_EHT, ITN=="no")$m_post_h_u.p)
boxplot(subset(Data_EHT, ITN=="LLIN")$m_post_h_u.p)
boxplot(subset(Data_EHT, ITN=="no")$f_h_u.p_live)
boxplot(subset(Data_EHT, ITN=="LLIN")$f_h_u.p_live)

hist(subset(Data_EHT, ITN=="no")$m_pre_h_u.p)
hist(subset(Data_EHT, ITN=="LLIN")$m_pre_h_u.p)
hist(subset(Data_EHT, ITN=="no")$m_post_h_u.p)
hist(subset(Data_EHT, ITN=="LLIN")$m_post_h_u.p)
hist(subset(Data_EHT, ITN=="no")$f_h_u.p_live)
hist(subset(Data_EHT, ITN=="LLIN")$f_h_u.p_live)


S = 0.9 # 0.7 to 0.95
g = 3 # 2,3,4,5
Pop = 1000 # 300 - 5000 (uniform)
prev = 0.5 # 0.1 - 0.8 (uniform)
k = 0.1 # 0.02 - 0.2 Churcher (but related to prevalence in human)
n = 11 # 8 - 15
HBI = 1 # 1

m_pre_h_u = 0.015 # 0 - 0.04 (EHT Moiroux)
m_pre_h_p = 0.5 # 0.1 - 1 (EHT Moiroux)           # TESTED
m_post_h_u = 0.005 # 0 - 0.01 (EHT Moiroux)
m_post_h_p = 0.21 # 0 - 0.5 (EHT Moiroux)
f_h_u_live = 0.55 # 0.2 - 0.8 (EHT Moiroux)
f_h_p_live = 0.31 # 0.1 - 0.7 (EHT Moiroux)       # TESTED
Ch = 0.6 # 0 - 1                                  # TESTED
pii = 0.9 # 0 - 1                                 # TESTED
pref_h_u = 0.5 # 0.3 - 0.9                        # TESTED