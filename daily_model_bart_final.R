##Assimilation of NEON observations into a process-based carbon cycle model reveals
##divergent mechanisms of carbon dynamics in temperate deciduous forests
##Chang Gyo Jung

##models start
##this model generates and updates covariance matrix for sampling
lapply(c("ggplot2","zoo","MASS","mvnfast","reshape2","lubridate","dplyr"), require, character.only = T)
rm(list = ls());graphics.off()
sapply(paste("sources/",list.files("sources/", pattern="*.R"),sep=""), source) #JC_JE; MCMC, ksi
input.data <- read.table("input/total_input.csv", header = T, sep = ",") #read input variables: PAR, RH, airT, soilT, VSWC
LMA <- read.table("input/LMA(BART_HARV_ORNL_STEI_UNDE).csv", header = T, sep = ",") #read LMA

site.list <- unique(input.data$siteID) #list of site

site <- site.list[1]

input.data <- input.data[input.data$siteID %in% site,]
LMA <- LMA[LMA$siteID %in% site,"LMA"]
attach(input.data) #data is available on the background
# detach(input.data) #data is available on the background

total.lengths <- length(airT)-365 #select 2016-2018
plot_no <- c(6,7) #

#VPD; from The ASCE Standardized Reference Evapotranspiration Equation
e_s <- 0.6108*exp((17.27*airT)/(airT+237.3))
e_a <- e_s*RH/100
VPD <- e_s-e_a
T_k <- airT+273.15 #C to K

GDD <- unlist(tapply(airT, year, GDD.func))
FDD <- unlist(tapply(airT, year, FDD.func))

#index
date <- Sys.Date()
No.test <- "8_2"
d <-15
##read observations
#NEE; leaf and wood litter; Fine root; SOC, microbial biomass, soil respiration (x)
observation <- read.table("input/observations(BART_HARV_ORNL_STEI_UNDE)_3.csv", header = T, sep =",")
obs.litter <- read.table("input/observations_litter(BART_HARV_ORNL_STEI_UNDE).csv", header = T, sep =",")

wood.init.2015 <- data.frame("siteID"=c("BART","HARV","ORNL","STEI"),"wood.2015"=c(14094.455,12280.721,22962.116,1510.171))
observation.w.o.2019 <- observation[observation$siteID %in% site&observation$year != 2019,] #selected years
obs.litter.w.o.2019 <- obs.litter[obs.litter$siteID %in% site & obs.litter$year.col != 2019,] #selected years

obs.SOC <- na.omit(observation.w.o.2019[,c("SOC","SOC.var","index")])
obs.NEE <- na.omit(observation.w.o.2019[,c("NEE.d", "index")])
obs.fine.root <- na.omit(observation.w.o.2019[,c("fine.root.mean","fine.root.var", "index")])
obs.leaf.litter <- obs.litter.w.o.2019[obs.litter.w.o.2019$functionalGroup %in% "Leaves",]
obs.wood.litter<- obs.litter.w.o.2019[obs.litter.w.o.2019$functionalGroup %in% "Woody material",]
obs.soilR<- na.omit(observation.w.o.2019[,c("soilR.d","soilR.d.var","index")])
obs.woody.biomass <- na.omit(observation.w.o.2019[,c("woody.biomass.mean","woody.biomass.var","index")])
obs.LAI <- na.omit(observation.w.o.2019[,c("LAI","LAI.var","index")])

obs.leaf.litter$litter.var <-
  ifelse(is.na(obs.leaf.litter$litter.var), obs.leaf.litter$litter.mean *0.8,obs.leaf.litter$litter.var)
obs.wood.litter$litter.var <-
  ifelse(is.na(obs.wood.litter$litter.var), obs.wood.litter$litter.mean *0.8,obs.wood.litter$litter.var)

###start!########

###MCMC
J_last <- 3000000000
pmin <- c(); pmax <- c() #create empty variables

####set prior ranges for
#GPP
#alpha_q  0.18  VM25  45        a1 6
pmin[1] <- 0.2; pmin[2] <- 40; pmin[3] <- 2
pmax[1] <- 0.5; pmax[2] <- 80; pmax[3] <- 20
#r_JmVm 1.7    E_k_c 129290       E_k_o 29395    Gamma_star_25 17
pmin[4] <- 1; pmin[5] <- 3.0e4; pmin[6] <-1.0e4;  pmin[7] <-10
pmax[4] <- 5; pmax[5] <- 1.5e5; pmax[6] <-6.0e4;  pmax[7] <-200
# E_gamma_25 66992  D_0  1.77     K_c_25    406       E_v_m  73310      K_o_25  0.32
pmin[8] <- 3.0e4; pmin[9] <- 1; pmin[10] <- 50; pmin[11] <- 1.0e4; pmin[12] <- 0.2
pmax[8] <- 1.0e5; pmax[9] <- 10; pmax[10] <- 600;pmax[11] <- 1.0e5; pmax[12] <- 0.5
#f_Ci (Hetheringtonand Woodward 2003))
##f_Ci 0.75     #f_Cs             #k_n 
pmin[13] <- 0.7;pmin[14] <- 0.80;pmin[15] <- 0.70 
pmax[13] <- 0.8;pmax[14] <- 0.95;pmax[15] <- 0.90

#GDD0   250     FDD0           
pmin[16] <-70; pmin[17] <- -4000;
pmax[16] <-100;pmax[17] <- -3000;

  #k1                   k2                    k3                    k4                  k5                k6                    k7
pmin[18] <- 5.00e-2; pmin[19] <- 3.20e-5; pmin[20] <- 9.20e-4; pmin[21] <- 8.57e-3; pmin[22] = 1.58e-3; pmin[23] = 2.50e-2; pmin[24] = 4.10e-5
pmax[18] <- 9.50e-2; pmax[19] <- 5.50e-5; pmax[20] <- 2.00e-3; pmax[21] <- 3.74e-2; pmax[22] = 4.50e-3; pmax[23] = 5.50e-2; pmax[24] = 1.24e-4

#offset.growing  
pmin[25] <- 0.4
pmax[25] <- 0.85

#matrix A 
#a41             #a43              #a64             a65            a75             a76              a67  
pmin[26] <- 0.5 ;pmin[27] <- 0.7; pmin[28] <- 0.4; pmin[29] = 0.1; pmin[30] = 0.2; pmin[31] = 0.1; pmin[32] = 0.1
pmax[26] <- 0.8 ;pmax[27] <- 0.9; pmax[28] <- 0.7; pmax[29] = 0.3; pmax[30] = 0.6; pmax[31] = 0.5; pmax[32] = 0.4

#matrix B
#b1                 b2               b3               
pmin[33] <- 0.1;pmin[34] <- 0.1;  pmin[35] <- 0.1;
pmax[33] <- 0.3;pmax[34] <- 0.3;  pmax[35] <- 0.3;

#Q10          
pmin[36] <- 1;
pmax[36] <- 4;

#initial pools
#x3root         x4f.lit         x5c.lit         x6mic.SOM       
pmin[37] <- 600;pmin[38] <- 100 ;pmin[39] <- 400;pmin[40] <- 80;
pmax[37] <- 700;pmax[38] <- 250;pmax[39] <- 500;pmax[40] <- 150;

#foliage.max  
pmin[41] <- 220
pmax[41] <- 280

########
pdif <- pmax-pmin; p_op <- (pmax+pmin)/2

#start simulation for generating covariance matrix
###initial pools
X <- matrix(0,7,total.lengths+1)
LAI_d <- c()
No.simu <- 100000

##for loop...
p_rec <- matrix(NA, length(pmin), No.simu)
p_upgraded <- matrix(NA, length(pmin), No.simu/2)
J_keep <- rep(NA,No.simu)

###first run#####
updated <- 0
sd <- 2.38/sqrt(length(pmin))

DJ <- list()
const.DJ <- 100000

#scaling
total.obs.no <- sum(length(obs.fine.root[,1]),
                    length(obs.leaf.litter[,1]),
                    length(obs.wood.litter[,1]),
                    length(obs.NEE[,1]),
                    length(obs.LAI[,1]),
                    length(obs.woody.biomass[,1]))

DJ.scale <- c((length(obs.fine.root[,1]))/total.obs.no,
              (length(obs.leaf.litter[,1]))/total.obs.no,
              (length(obs.wood.litter[,1]))/total.obs.no,
              (length(obs.NEE[,1]))/total.obs.no,
              (length(obs.LAI[,1]))/total.obs.no,
              (length(obs.woody.biomass[,1]))/total.obs.no)

DJ[1] <- list(2*obs.fine.root$fine.root.var) #root
DJ[2] <- list(2*obs.leaf.litter$litter.var*100) #leaf.litter
DJ[3] <- list(2*obs.wood.litter$litter.var*20) #wood.litter
DJ[4] <- list(2*(0.04*obs.NEE$NEE.d)^2) #NEE
DJ[5] <- list(2*obs.LAI$LAI.var*0.0001) #LAI
DJ[6] <- list(2*obs.woody.biomass$woody.biomass.var) #woody biomass

for (i in 1:length(DJ)) {
  DJ[[i]] <- DJ[[i]]*DJ.scale[i]
}

for (simu in 1:No.simu) {
  ##update pars
  if (simu <= 10000) {
    p_new <- generate.pars(p_op, pmin, pmax, d)
  } else {
      p_new <- generate.pars.cov(p_op, covars)
      }
  
  p_rec[,simu] <- p_new
  par.names <- c("alpha_q","VM25","a1","r_JmVm","E_k_c","E_k_o","Gamma_star_25",
                 "E_gamma_25","D_0","K_c_25","E_v_m","K_o_25",
                 "f_Ci","f_Cs","k_n","GDD0","FDD0",           
                 paste("k",1:7,sep=""),"offset.growing","a41", "a43","a64","a65","a75","a76","a67","b1","b2","b3","Q10",
                 paste("X",3:6,sep = ""),"foliage.max")
  for (i in 1:length(par.names)) {
    assign(par.names[i], p_new[i])
  }
  
  #Photosynthesis
  o_x <- 0.21 #;           %mol(O2)/mol(air)).
  C_a <- 410 #;            % air CO2 concentration
  
  #assign pars for:
  k_comp <- c(k1,k2,k3,k4,k5,k6,k7) #for matrix K 
  a_comp <- c(-1,      0,       0,     0,    0,  0,   0,
               0,     -1,       0,     0,    0,  0,   0,
               0,      0,      -1,     0,    0,  0,   0,
              a41,      0,     a43,    -1,    0,  0,   0,
            1-a41,      1,   1-a43,     0,   -1,  0,   0,
               0,      0,       0,   a64,  a65, -1, a67,
               0,      0,       0,     0,  a75, a76, -1)

  K_mat <- diag(k_comp,7,7)
  A_mat <- matrix(a_comp, nrow = 7, ncol = 7, byrow = T)
  B_mat <- rbind(b1, b2, b3, 0, 0, 0, 0) #vector B 
  
  #environmental scalar
  t.m.adj <- ksi(total.lengths, soilT, VSWC, 0.2,Q10) 
  
  LAI <- 0.1 #the first day "LAI"  
  LAI_d[1] <- 0.1 #the first day
  
  #initial pools
  X2 <- wood.init.2015[wood.init.2015$siteID %in% site,"wood.2015"]
  X7 <- 9273.569 - X6
  x0 <- cbind(c(0,X2,X3,X4,X5,X6,X7))
  X[,1] <- x0
  
  #switch; phenology 
  leaf.sw <- ifelse((GDD < GDD0) | (FDD < FDD0), 1,offset.growing) 
  
  C_i <- C_a * f_Ci 
  C_s <- C_a * f_Cs #surface co2
  
  #simulation An
  GPP_model <- c()
  Gamma_star <- Gamma_star_25*exp(E_gamma_25*(T_k-298)/(8.314*298.15*T_k))
  Vcmax <- VM25*exp(E_v_m*(T_k-298)/(8.314*298.15*T_k))
  J_c <- JC_1(Gamma_star,T_k,VM25,E_v_m,K_c_25,E_k_c,K_o_25,E_k_o,o_x,Vcmax,C_i)
  J_e <- JE_1(Gamma_star,alpha_q,PAR,C_i,r_JmVm,VM25,T_k,E_v_m)
  A <- pmin(J_c, J_e) #Farquhar 
  A <- pmax(A, 0.01)
  Gs <- 0.01 + a1*A/((C_s-Gamma_star)*(1+VPD/D_0)) #Leuning stomatal conductance
  An <- Gs*(C_a - C_i) #carbon influx of the top leaf layer
  
  for (i in 1:total.lengths) {
    Ac <- An[i]*(1-exp(-k_n*LAI))/k_n #scale-up to canopy
    GPP_model[i] <- Ac*1.037750 #convert g C m-2 d-1 
    t.m.adj.mat <- diag(c(1,1,1,rep(t.m.adj[i],4))) #environmental scalar 
    K_mat.2 <- K_mat %*% diag(c(leaf.sw[i],rep(1,6)),7,7) #matrix K with switch
    max.leaf.sw <- ifelse(X[1,i]>foliage.max, 0, 1)
    
    #update X
    X[,i+1] <- (diag(1,7,7) + A_mat %*% K_mat.2 %*% t.m.adj.mat) %*% X[,i] + diag(c(max.leaf.sw,rep(1,6)),7,7)%*%B_mat%*%GPP_model[i]
    LAI <- max((X[1,i+1]/LMA),0.1)
    LAI_d[i+1] <- LAI
  }  

  X.2 <- X[,1:total.lengths] 
  
  #mapping function
  phi_leaf <- c(1, 0, 0, 0, 0, 0, 0)
  phi_wood <- c(0, 1, 0, 0, 0, 0, 0)
  phi_het_Resp <- c(0,0,0,(1-a64)*k4, (1-a65-a75)*k5, (1-a76)*k6, (1-a67)*k7)
  phi_root <- c(0, 0, 1, 0, 0, 0, 0)
  phi_fine_litter  <- c(1*k1, 0, 0, 0, 0, 0, 0)
  phi_coarse_litter <- c(0, 1*k2, 0, 0,0, 0, 0)
  phi_SOM  = c(0, 0, 0, 0, 0, 1, 1)
  
  #simulation for comparing to observation 
  simu.leaf <- phi_leaf %*% X.2
  simu.wood <- phi_wood %*% X.2
  simu.root <- phi_root %*% X.2
  simu.leaf.litter.all <- phi_fine_litter %*% X.2
  simu.wood.litter.all <- phi_coarse_litter %*% X.2
  
  #for litter
  litter.index.b.leaf <- obs.leaf.litter$index.set
  litter.index.b.wood <- obs.wood.litter$index.set
  litter.index.leaf <- obs.leaf.litter$index.col
  litter.index.wood <- obs.wood.litter$index.col
  simu.leaf.litter <- c()
  simu.wood.litter <- c()
  
  for (i in 1:length(litter.index.b.leaf)) {
    simu.leaf.litter[i] <- sum(simu.leaf.litter.all[litter.index.b.leaf[i]:litter.index.leaf[i]])
  }

  for (i in 1:length(litter.index.b.wood)) {
    simu.wood.litter[i] <- sum(simu.wood.litter.all[litter.index.b.wood[i]:litter.index.wood[i]])
  }
  
  simu.leaf.litter <- as.matrix(simu.leaf.litter)
  simu.wood.litter <- as.matrix(simu.wood.litter)
  simu.soc <- phi_SOM %*% X.2
  simu.autoR <- (1-sum(B_mat[1:3]))*GPP_model
  simu.hetR <- t.m.adj*(phi_het_Resp %*% X.2)[1,]
  
  simu.soilR <-  simu.hetR + simu.autoR*(B_mat[3]/sum(B_mat))
  simu.NEE <- as.matrix(simu.autoR + simu.hetR - GPP_model)

  #calculate J
  J <- list()
  J_new <- c()
  J[1] <- list((simu.root[,obs.fine.root$index] - obs.fine.root$fine.root.mean)^2)
  J[2] <- list((simu.leaf.litter - obs.leaf.litter$litter.mean)^2) # the correct way to calculate errors in the cost function
  J[3] <- list((simu.wood.litter - obs.wood.litter$litter.mean)^2)
  J[4] <- list((simu.NEE[obs.NEE$index] - obs.NEE$NEE.d)^2)
  J[5] <- list((LAI_d[obs.LAI$index] - obs.LAI$LAI)^2)
  J[6] <- list((simu.wood[obs.woody.biomass$index] - obs.woody.biomass$woody.biomass.mean)^2)
  
  for (i in 1:length(J)) {
    J_new[i] <- sum(J[[i]]/(DJ[[i]]*const.DJ), na.rm=T)
  }
  J_new2 <- sum(J_new, na.rm = T) 
  delta.J <- J_new2 - J_last
  # print(J_new)
  if (min(1, exp(-delta.J)) > runif(1)) {
    p_op <- p_new
    J_last <- J_new2
    updated <- updated + 1
    J_keep[updated] <- J_new2 
    p_upgraded[,updated] <- p_op
    
  }
  print(paste("simu =", simu, "updated =", updated))
  
  if (simu == 10000) {
    covars <- cov(t(p_rec[,1:simu]))
  }
  if (simu > 10000) {
    covars <- sd*cov(t(p_rec[,1:simu]))
  }
}

par(mfrow=plot_no)
par(mar=c(1.8,1.7,1.5,1.5))

for (i in 1:length(p_upgraded[,1])) {
  hist(p_upgraded[i,1:updated], xlab =NA, ylab = NA, main = par.names[i], xlim=c(pmin[i],pmax[i]))
}

par(mfrow=plot_no)
par(mar=c(1.8,1.8,1,1))
for (i in 1:length(p_upgraded[,1])) {
  plot(p_upgraded[i,1:updated],type="l", xlab =NA, ylab = NA, main = par.names[i])
}

png(paste("pars_all", date, site, No.test, ".png", sep =  "_"), width = 1000, height = 1000, units = "px")
par(mfrow=plot_no)
par(mar=c(1.8,1.8,1,1))
for (i in 1:length(p_upgraded[,1])) {
  plot(p_upgraded[i,round(updated/2):updated],type="l", xlab =NA, ylab = NA, main = par.names[i])
}
dev.off()

png(paste("pars_hist", date, site,  No.test, ".png", sep =  "_"), width = 1000, height = 1000, units = "px")
par(mfrow=plot_no)
par(mar=c(1.8,1.8,1,1))
for (i in 1:length(p_upgraded[,1])) {
  hist(p_upgraded[i,round(updated/2):updated], xlab =NA, ylab = NA, main = par.names[i], xlim=c(pmin[i],pmax[i]))
}
dev.off()

p_upgraded.plot <- reshape2::melt(t(p_upgraded[,round(updated/2):updated]))
p_upgraded.plot$Var2 <- factor(p_upgraded.plot$Var2, labels = par.names)
hist_par <- ggplot(p_upgraded.plot, aes(value)) +
  geom_histogram(bins = 50) +
  facet_wrap(.~Var2, scales = "free") +
  theme_bw() +
  scale_x_continuous(limits = c()) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1),
        panel.grid = element_blank())
ggsave(paste("hist_par", date, site, No.test,".png",sep = "_"), hist_par, width = 22, height = 22, units = "cm", device = "png")

write.csv(p_upgraded[,round(updated/2):updated],paste("pars_updated",date, site, No.test, ".csv",sep="_"), row.names = F)

sink(paste(date,site, No.test, "pars.txt", sep = "_"));print(rbind(par.names, pmin, pmax));sink()

# #######
# ##simul
# #######
# # 200 sample
LAI.ac <- c()
simu.wood <- c()
simu.root <- c()
simu.leaf.litter.ac <- c()
simu.wood.litter.ac <- c()
simu.soc <- c()
simu.autoR <- c()
simu.hetR <- c()
simu.er <- c()
simu.soilR <-  c()
simu.NEE <- c()
GPP_model.ac <- c()
GDD.all <- c()
FDD.all <- c()

simu.leaf.p.ac <- c()
simu.leaf.litter.p.ac <- c()
simu.wood.litter.p.ac <- c()
simu.fast.p.ac <- c()
simu.slow.p.ac <- c()

#load data again
total.lengths <- length(airT)
plot_no <- c(6,7)

observation <- read.table("input/observations(BART_HARV_ORNL_STEI_UNDE)_3.csv", header = T, sep =",")
obs.litter <- read.table("input/observations_litter(BART_HARV_ORNL_STEI_UNDE).csv", header = T, sep =",")

observation <- observation[observation$siteID %in% site,]
obs.litter <- obs.litter[obs.litter$siteID %in% site,] 

##read observations
obs.SOC <- na.omit(observation[,c("SOC","SOC.var", "index")])
obs.NEE <- na.omit(observation[,c("NEE.d", "index")])
obs.fine.root <- na.omit(observation[,c("fine.root.mean","fine.root.var", "index")])
obs.leaf.litter <- obs.litter[obs.litter$functionalGroup %in% "Leaves",]
obs.wood.litter<- obs.litter[obs.litter$functionalGroup %in% "Woody material",]
obs.soilR<- na.omit(observation[,c("soilR.d","soilR.d.var","index")])
obs.woody.biomass <- na.omit(observation[,c("woody.biomass.mean","woody.biomass.var","index")])
obs.LAI <- na.omit(observation[,c("LAI","LAI.var","index")])

X <- matrix(0,7,total.lengths+1)

sample.index <- sample(round(updated/2):updated,200, replace = T)
p_new2 <- c()
p_new2 <- p_upgraded[,sample.index]

data.mle <- c()
for (i in 1:length(p_upgraded[,1])) {
  data.mle[i] <-  MASS::fitdistr(p_upgraded[i,round(updated/2):updated], "normal")$estimate[[1]]
}

p_new2 <- cbind(p_new2,data.mle)

for (j in 1:200) {
  p_new <- p_new2[,j]
  for (i in 1:length(par.names)) {
    assign(par.names[i], p_new[i])
  }
  #assign pars for:
  k_comp <- c(k1,k2,k3,k4,k5,k6,k7) #matrix K
  a_comp <- c(-1,     0,       0,     0,    0,  0,   0,
              0,     -1,       0,     0,    0,  0,   0,
              0,      0,      -1,     0,    0,  0,   0,
              a41,    0,     a43,    -1,    0,  0,   0,
              1-a41, 1,   1-a43,     0,   -1,  0,   0,
              0,      0,       0,   a64,  a65, -1, a67,
              0,      0,       0,     0,  a75, a76, -1)
  
  K_mat <- diag(k_comp,7,7)
  A_mat <- matrix(a_comp, nrow = 7, ncol = 7, byrow = T)
  B_mat <- rbind(b1, b2, b3, 0, 0, 0, 0) #matrix B
  #environmental scalar
  t.m.adj <- ksi(total.lengths, soilT, VSWC, 0.2,Q10) 
  
  LAI <- 0.1 #the first day "LAI"  
  LAI_d[1] <- 0.1 #the first day
  
  #initial pools
  X2 <- wood.init.2015[wood.init.2015$siteID %in% site,"wood.2015"]
  x0 <- cbind(c(0,X2,X3,X4,X5,X6,X7))
  X[,1] <- x0
  #switch; phenology 
  leaf.sw <- ifelse((GDD < GDD0) | (FDD < FDD0), 1,offset.growing) 
  
#Photosynthesis
C_i <- C_a * f_Ci 
C_s <- C_a * f_Cs #surface co2

#simulation An
GPP_model <- c()
Gamma_star <- Gamma_star_25*exp(E_gamma_25*(T_k-298)/(8.314*298.15*T_k))
Vcmax <- VM25*exp(E_v_m*(T_k-298)/(8.314*298.15*T_k))
J_c <- JC_1(Gamma_star,T_k,VM25,E_v_m,K_c_25,E_k_c,K_o_25,E_k_o,o_x,Vcmax,C_i)
J_e <- JE_1(Gamma_star,alpha_q,PAR,C_i,r_JmVm,VM25,T_k,E_v_m)
A <- pmin(J_c, J_e) #Farquhar 
A <- pmax(A, 0.01)
Gs <- 0.01 + a1*A/((C_s-Gamma_star)*(1+VPD/D_0)) #Leuning stomatal conductance
An <- Gs*(C_a - C_i) #carbon influx of the top leaf layer
t.m.adj <- ksi(total.lengths, soilT, VSWC, 0.2,Q10) 

LAI[1] <- 0.1
LAI_d[1] <- 0.1
#initial pools
X2 <- wood.init.2015[wood.init.2015$siteID %in% site,"wood.2015"]
X7 <- 9273.569 - X6 
x0 <- cbind(c(0,X2,X3,X4,X5,X6,X7))
X[,1] <- x0

for (i in 1:total.lengths) {
  Ac <- An[i]*(1-exp(-k_n*LAI))/k_n #scale-up to canopy
  GPP_model[i] <- Ac*1.037750 #convert g C m-2 d-1 
  t.m.adj.mat <- diag(c(1,1,1,rep(t.m.adj[i],4))) #environmental scalar 
  K_mat.2 <- K_mat %*% diag(c(leaf.sw[i],rep(1,6)),7,7) #matrix K with switch
  max.leaf.sw <- ifelse(X[1,i]>foliage.max, 0, 1)
  #update X
  X[,i+1] <- (diag(1,7,7) + A_mat %*% K_mat.2 %*% t.m.adj.mat) %*% X[,i] + diag(c(max.leaf.sw,rep(1,6)),7,7)%*%B_mat%*%GPP_model[i]
  LAI <- max((X[1,i+1]/LMA),0.1)
  LAI_d[i+1] <- LAI
}  

X.2 <- X[,1:total.lengths] 

#mapping function
phi_leaf <- c(1, 0, 0, 0, 0, 0, 0)
phi_wood <- c(0, 1, 0, 0, 0, 0, 0)
phi_het_Resp <- c(0,0,0,(1-a64)*k4, (1-a65-a75)*k5, (1-a76)*k6, (1-a67)*k7)
phi_root <- c(0, 0, 1, 0, 0, 0, 0)
phi_fine_litter  <- c(1*k1, 0, 0, 0, 0, 0, 0)
phi_coarse_litter <- c(0, 1*k2, 0, 0,0, 0, 0)
phi_SOM  = c(0, 0, 0, 0, 0, 1, 1)

#for litter
simu.leaf.litter.all <- phi_fine_litter %*% X.2
simu.wood.litter.all <- phi_coarse_litter %*% X.2
litter.index.b.leaf <- obs.leaf.litter$index.set
litter.index.b.wood <- obs.wood.litter$index.set
litter.index.leaf <- obs.leaf.litter$index.col
litter.index.wood <- obs.wood.litter$index.col
simu.leaf.litter <- c()
simu.wood.litter <- c()

for (i in 1:length(litter.index.b.leaf)) {
  simu.leaf.litter[i] <- sum(simu.leaf.litter.all[litter.index.b.leaf[i]:litter.index.leaf[i]])
}

for (i in 1:length(litter.index.b.wood)) {
  simu.wood.litter[i] <- sum(simu.wood.litter.all[litter.index.b.wood[i]:litter.index.wood[i]])
}

simu.leaf.litter <- as.matrix(simu.leaf.litter)
simu.wood.litter <- as.matrix(simu.wood.litter)

###gethering
LAI.ac[[j]] <- LAI_d
GPP_model.ac[[j]] <- GPP_model
simu.wood[[j]] <- phi_wood %*% X.2
simu.root[[j]] <- phi_root %*% X.2
simu.leaf.litter.ac[[j]] <- as.matrix(simu.leaf.litter)
simu.wood.litter.ac[[j]] <- as.matrix(simu.wood.litter)
simu.soc[[j]] <- phi_SOM %*% X.2
simu.autoR[[j]] <- (1-sum(B_mat[1:3]))*GPP_model
simu.hetR[[j]] <- t.m.adj*(phi_het_Resp %*% X.2)[1,]
simu.er[[j]] <- simu.autoR[[j]] + simu.hetR[[j]]
simu.soilR[[j]] <- t.m.adj*(phi_het_Resp %*% X.2)[1,] + simu.autoR[[j]]*(B_mat[3]/sum(B_mat))
simu.NEE[[j]] <- as.matrix(simu.autoR[[j]] + simu.hetR[[j]] - GPP_model)

simu.leaf.p.ac[[j]] <- X.2[1,]
simu.leaf.litter.p.ac[[j]] <- X.2[4,]
simu.wood.litter.p.ac[[j]] <- X.2[5,]
simu.fast.p.ac[[j]] <- X.2[6,]
simu.slow.p.ac[[j]] <- X.2[7,]

#calculate onset, offset dates
GDDD <- ifelse((GDD < GDD0), 1,NA)
FDDD <- ifelse((FDD < FDD0), NA,1) 
GDD.all[[j]] <- as.numeric(tapply(GDDD, year, sum, na.rm = T)+1)
FDD.all[[j]] <- as.numeric(tapply(FDDD, year, sum, na.rm = T)+1)
}

GDD.mean <- round(colMeans(do.call(rbind,GDD.all)))
FDD.mean <- round(colMeans(do.call(rbind,FDD.all)))

LAI.all <- t(do.call(rbind,LAI.ac))[1:total.lengths,]
simu.wood.all <- t(do.call(rbind,simu.wood))
simu.root.all <- t(do.call(rbind,simu.root))
simu.soc.all <- t(do.call(rbind,simu.soc))
simu.NEE.all <- (do.call(cbind,simu.NEE))
simu.GPP.all <- (do.call(cbind,GPP_model.ac))
simu.soilR.all <- (do.call(cbind,simu.soilR))
simu.autoR.all <- (do.call(cbind,simu.autoR))
simu.hetR.all <- (do.call(cbind,simu.hetR))
simu.ER.all <- (do.call(cbind,simu.er))
simu.leaf.litter.ac.all <- do.call(cbind,simu.leaf.litter.ac)
simu.wood.litter.ac.all <- do.call(cbind,simu.wood.litter.ac)

#save
save.image(paste(date,site,No.test, "_2.RData", sep = "_"))
