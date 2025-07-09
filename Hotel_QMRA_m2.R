
#Model 2: Susceptible touch surfaces 4hour after an infected guest passes by the hotel lobby


fomite <-"Elevator"
source("Hotel_Parameters.R")


##M2.S1: "Elevator" touch____________________________________________   

###M2.S1.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.21<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.21)<-eventsname

Conc.s.21<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.21)<-eventsname

Dose.21<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.21)<-eventsname

Risk.21<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.21)<-eventsname

#Event 0. Infected person passed by hotel lobby

##Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
fomite <-"Elevator"
source("Hotel_Parameters.R")
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf


#Event 1. Susceptible guest contact seq 1 (Elevator)

##Concentration of norovirus on (Elevator) after 4 hours of ordinary hotel lobby usuage
Conc.4hour<-(Conc.recover/(Conc.seed*RE.swab))*Conc.surf.i

Conc.h.21[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.4hour*exp(-k.surf*Time.m1))}
Conc.s.21[1,]<-Conc.4hour*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.4hour*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.21[1,]<-0
Risk.21[1,]<-0

#Event 2. Hand to face contact
Conc.h.21[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.21[1,]*exp(-k.hand*Time.m1))
Conc.s.21[2,]<-Conc.s.21[1,]*exp(-k.surf*Time.m1)
Dose.21[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.21[1,]*exp(-k.hand*Time.m1)
Risk.21[2,]<-safe_risk_calc(Dose.21[2,])


##M2.S2: "Front desk" touch__________________________________


###M2.S2.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.22<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.22)<-eventsname

Conc.s.22<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.22)<-eventsname

Dose.22<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.22)<-eventsname

Risk.22<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.22)<-eventsname

#Event 0. Infected person passed by hotel lobby

##Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
fomite <-"Elevator"
source("Hotel_Parameters.R")
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf


#Event 1. Susceptible guest contact seq 1 (Front desk)

##Concentration of norovirus on (Front desk) after 4 hours of ordinary hotel lobby usuage
fomite <-"Frontdesk"
source("Hotel_Parameters.R")
Conc.4hour<-(Conc.recover/(Conc.seed*RE.swab))*Conc.surf.i

Conc.h.22[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.4hour*exp(-k.surf*Time.m1))}
Conc.s.22[1,]<-Conc.4hour*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.4hour*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.22[1,]<-0
Risk.22[1,]<-0

#Event 2. Hand to face contact
Conc.h.22[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.22[1,]*exp(-k.hand*Time.m1))
Conc.s.22[2,]<-Conc.s.22[1,]*exp(-k.surf*Time.m1)
Dose.22[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.22[1,]*exp(-k.hand*Time.m1)
Risk.22[2,]<-safe_risk_calc(Dose.22[2,])


##M2.S3: "Elevator - Front desk" touch__________________________________
###M2.S3.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.23<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.23)<-eventsname

Conc.s.23<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.23)<-eventsname

Dose.23<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.23)<-eventsname

Risk.23<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.23)<-eventsname

#Event 0. Infected person passed by hotel lobby

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
fomite <-"Elevator"
source("Hotel_Parameters.R")
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf

#Concentration of norovirus on surface i after 4 hours of ordinary hotel lobby usuage
Conc.4hour<-(Conc.recover/(Conc.seed*RE.swab))*Conc.surf.i

#Event 1. Susceptible guest contact seq 1 (Elevator)
Conc.h.23[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.4hour*exp(-k.surf*Time.m1))}
Conc.s.23[1,]<-Conc.4hour*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.4hour*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.23[1,]<-0
Risk.23[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Front desk)
fomite <-"Frontdesk"
source("Hotel_Parameters.R")
#Concentration of norovirus on surface i after 4 hours of ordinary hotel lobby usuage
Conc.4hour<-(Conc.recover/(Conc.seed*RE.swab))*Conc.surf.i

Conc.h.23[2,]<-Conc.h.23[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.23[1,]*exp(-k.hand*Time.m1)-TE.sh*Conc.4hour*exp(-k.surf*2*Time.m1))}
Conc.s.23[2,]<-Conc.4hour*exp(-k.surf*2*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.4hour*exp(-k.surf*2*Time.m1)-Conc.h.23[1,]*exp(-k.hand*Time.m1))}
Dose.23[2,]<-0
Risk.23[2,]<-0


#Event 3. Hand to face contact
Conc.h.23[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.23[2,]*exp(-k.hand*Time.m1))
Conc.s.23[3,]<-Conc.s.23[2,]*exp(-k.surf*Time.m1)
Dose.23[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.23[2,]*exp(-k.hand*Time.m1)
Risk.23[3,]<-safe_risk_calc(Dose.23[3,])


##M2.S4: "Table top" touch__________________________________

###M2.S4.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.24<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.24)<-eventsname

Conc.s.24<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.24)<-eventsname

Dose.24<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.24)<-eventsname

Risk.24<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.24)<-eventsname

#Event 0. Infected person passed by hotel lobby

##Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
fomite <-"Elevator"
source("Hotel_Parameters.R")
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf


#Event 1. Susceptible guest contact seq 1 (Table top)

##Concentration of norovirus on (Table top) after 4 hours of ordinary hotel lobby usuage
fomite <-"Table"
source("Hotel_Parameters.R")
Conc.4hour<-(Conc.recover/(Conc.seed*RE.swab))*Conc.surf.i

Conc.h.24[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.4hour*exp(-k.surf*Time.m1))}
Conc.s.24[1,]<-Conc.4hour*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.4hour*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.24[1,]<-0
Risk.24[1,]<-0

#Event 2. Hand to face contact
Conc.h.24[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.24[1,]*exp(-k.hand*Time.m1))
Conc.s.24[2,]<-Conc.s.24[1,]*exp(-k.surf*Time.m1)
Dose.24[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.24[1,]*exp(-k.hand*Time.m1)
Risk.24[2,]<-safe_risk_calc(Dose.24[2,])


##M2.S5: "Elevator - Table top" touch___________________________________

###M2.s5.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.25<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.25)<-eventsname

Conc.s.25<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.25)<-eventsname

Dose.25<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.25)<-eventsname

Risk.25<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.25)<-eventsname

#Event 0. Infected person passed by hotel lobby
fomite <-"Elevator"
source("Hotel_Parameters.R")

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf


#Event 1. Susceptible guest contact seq 1 (Elevator)
#Concentration of norovirus on surface i after 4 hours of ordinary hotel lobby usuage
Conc.4hour<-(Conc.recover/(Conc.seed*RE.swab))*Conc.surf.i

Conc.h.25[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.4hour*exp(-k.surf*Time.m1))}
Conc.s.25[1,]<-Conc.4hour*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.4hour*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.25[1,]<-0
Risk.25[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Table))
fomite <-"Table"
source("Hotel_Parameters.R")
#Concentration of norovirus on surface i after 4 hours of ordinary hotel lobby usuage
Conc.4hour<-(Conc.recover/(Conc.seed*RE.swab))*Conc.surf.i


Conc.h.25[2,]<-Conc.h.25[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.25[1,]*exp(-k.hand*Time.m1)-TE.sh*Conc.4hour*exp(-k.surf*2*Time.m1))}
Conc.s.25[2,]<-Conc.4hour*exp(-k.surf*2*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.4hour*exp(-k.surf*2*Time.m1)-Conc.h.25[1,]*exp(-k.hand*Time.m1))}
Dose.25[2,]<-0
Risk.25[2,]<-0


#Event 3. Hand to face contact
Conc.h.25[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.25[2,]*exp(-k.hand*Time.m1))
Conc.s.25[3,]<-Conc.s.25[2,]*exp(-k.surf*Time.m1)
Dose.25[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.25[2,]*exp(-k.hand*Time.m1)
Risk.25[3,]<-safe_risk_calc(Dose.25[3,])


##5. plotting------------------------------------------------------------
library(ggplot2)
library(ggpubr)


Risk.21.df<-as.data.frame(t(Risk.21))
Risk.22.df<-as.data.frame(t(Risk.22))
Risk.23.df<-as.data.frame(t(Risk.23))
Risk.24.df<-as.data.frame(t(Risk.24))
Risk.25.df<-as.data.frame(t(Risk.25))

type<-c(rep("E",iter),rep("FD",iter),rep("E-FD",iter),rep("TT",iter),rep("E-TT",iter))
value<-c(Risk.21.df$"Hand to face touch",Risk.22.df$"Hand to face touch", Risk.23.df$"Hand to face touch", Risk.24.df$"Hand to face touch", Risk.25.df$"Hand to face touch")  

data.m2<-data.frame(type,value)

windows()
p<- ggplot(data.m2)+geom_violin(aes(x=type,y=value,fill=type, group=type),alpha=0.3,draw_quantiles = c(0.25,0.5,0.75))
p+ scale_y_continuous(trans="log10", breaks=c(1e-8, 1e-10, 1e-12, 1e-14, 1e-16, 1e-18, 1e-20),
                      label=c(-8, -10, -12, -14, -16, -18, -20)) +
  scale_x_discrete(limits=c("E","FD","E-FD","TT","E-TT"))+
  ggtitle("Model 2: Susceptible guest touches fomite 4 hours after lobby use by infected guest")+
  xlab("Susceptible Guest Touch Sequences")+
  ylab(expression("Infection risk of norovirus (log" [10] * ")"))

ggsave("Model2.tiff", dpi=600, dev='tiff', height=6, width=9, units="in")

### Debugging ---------------------------
# 각 fomite별 Conc.4hour 값 확인
fomite_types <- c("Elevator", "Frontdesk", "Table")
for(ft in fomite_types) {
  fomite <- ft
  source("Hotel_Parameters.R")
  
  Conc.surf.i <- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf
  Conc.4hour <- (Conc.recover/(Conc.seed*RE.swab))*Conc.surf.i
  
  cat(ft, ":\n")
  cat("  Conc.recover:", range(Conc.recover), "\n")
  cat("  Conc.4hour 범위:", range(Conc.4hour), "\n")
  cat("  0인 값 개수:", sum(Conc.4hour == 0), "/", length(Conc.4hour), "\n")
  cat("  매우 작은 값(<1e-50) 개수:", sum(Conc.4hour < 1e-50 & Conc.4hour > 0), "\n\n")
}
# Risk 계산 과정을 단계별로 확인
fomite_types <- c("Frontdesk", "Table")
scenarios <- list("Risk.22" = 22, "Risk.24" = 24)

for(i in 1:length(fomite_types)) {
  ft <- fomite_types[i]
  scenario_num <- scenarios[[i]]
  
  cat("=== ", ft, " (Risk.", scenario_num, ") ===\n")
  
  # 초기 설정
  fomite <- "Elevator"
  source("Hotel_Parameters.R")
  Conc.surf.i <- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf
  
  # 해당 fomite로 변경
  fomite <- ft
  source("Hotel_Parameters.R")
  Conc.4hour <- (Conc.recover/(Conc.seed*RE.swab))*Conc.surf.i
  
  # 1단계: 손 농도 계산
  if(ft == "Frontdesk") {
    Conc.h.step1 <- Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.4hour*exp(-k.surf*Time.m1))}
  } else { # Table
    Conc.h.step1 <- TE.sh*Conc.4hour*exp(-k.surf*Time.m1)
  }
  
  # 2단계: 최종 손 농도 (hand-to-face 직전)
  Conc.h.final <- Conc.h.step1*exp(-k.hand*Time.m1)
  
  # 3단계: Dose 계산
  Dose.final <- TE.hf*Frac.hf*T.handarea*Conc.h.final
  
  # 4단계: Risk 계산
  Risk.final <- safe_risk_calc(Dose.final)
  
  cat("  Conc.h.step1 - 0인 값:", sum(Conc.h.step1 == 0), "/", length(Conc.h.step1), "\n")
  cat("  Conc.h.final - 0인 값:", sum(Conc.h.final == 0), "/", length(Conc.h.final), "\n")
  cat("  Dose.final - 0인 값:", sum(Dose.final == 0), "/", length(Dose.final), "\n")
  cat("  Risk.final - 0인 값:", sum(Risk.final == 0), "/", length(Risk.final), "\n")
  cat("  Conc.h.step1 범위:", range(Conc.h.step1), "\n")
  cat("  Dose.final 범위:", range(Dose.final), "\n\n")
}

##6. Data pulling--(Need to---------------------------------------------------------
#
##6.1 2ACH data 

#Conc
matrix.Conc_2<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Conc_2)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Conc_2)<-studentsname


for (f in 1:numstudents){
  matrix.Conc_2[f,1]<-mean(Conc_2[f,])
  matrix.Conc_2[f,2]<-sd(Conc_2[f,])
  matrix.Conc_2[f,3]<-min(Conc_2[f,])
  matrix.Conc_2[f,4]<-max(Conc_2[f,])
}

#Dose
matrix.Dose_2<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Dose_2)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Dose_2)<-studentsname

for (f in 1:numstudents){
  matrix.Dose_2[f,1]<-mean(Dose_2[f,])
  matrix.Dose_2[f,2]<-sd(Dose_2[f,])
  matrix.Dose_2[f,3]<-min(Dose_2[f,])
  matrix.Dose_2[f,4]<-max(Dose_2[f,])
}

#Risk 
matrix.Risk_2<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Risk_2)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Risk_2)<-studentsname

for (f in 1:numstudents){
  matrix.Risk_2[f,1]<-mean(Risk_2[f,])
  matrix.Risk_2[f,2]<-sd(Risk_2[f,])
  matrix.Risk_2[f,3]<-min(Risk_2[f,])
  matrix.Risk_2[f,4]<-max(Risk_2[f,])
}
#Check the data 
View(matrix.Conc_2)
View(matrix.Dose_2)
View(matrix.Risk_2)

#Pull out the data 
library(openxlsx)
write.csv(matrix.Conc_2, file="Conc_2.csv")
write.csv(matrix.Dose_2, file="Dose_2.csv")
write.csv(matrix.Risk_2, file="Risk_2.csv")

##6.2 4ACH data------------------------------------------------------------

#Conc
matrix.Conc_4<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Conc_4)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Conc_4)<-studentsname


for (f in 1:numstudents){
  matrix.Conc_4[f,1]<-mean(Conc_4[f,])
  matrix.Conc_4[f,2]<-sd(Conc_4[f,])
  matrix.Conc_4[f,3]<-min(Conc_4[f,])
  matrix.Conc_4[f,4]<-max(Conc_4[f,])
}

#Dose
matrix.Dose_4<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Dose_4)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Dose_4)<-studentsname

for (f in 1:numstudents){
  matrix.Dose_4[f,1]<-mean(Dose_4[f,])
  matrix.Dose_4[f,2]<-sd(Dose_4[f,])
  matrix.Dose_4[f,3]<-min(Dose_4[f,])
  matrix.Dose_4[f,4]<-max(Dose_4[f,])
}

#Risk 
matrix.Risk_4<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Risk_4)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Risk_4)<-studentsname

for (f in 1:numstudents){
  matrix.Risk_4[f,1]<-mean(Risk_4[f,])
  matrix.Risk_4[f,2]<-sd(Risk_4[f,])
  matrix.Risk_4[f,3]<-min(Risk_4[f,])
  matrix.Risk_4[f,4]<-max(Risk_4[f,])
}
#Check the data 
View(matrix.Conc_4)
View(matrix.Dose_4)
View(matrix.Risk_4)

#Pull out the data 
library(openxlsx)
write.csv(matrix.Conc_4, file="Conc_4.csv")
write.csv(matrix.Dose_4, file="Dose_4.csv")
write.csv(matrix.Risk_4, file="Risk_4.csv")

##6.3 6ACH data------------------------------------------- 

#Conc
matrix.Conc_6<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Conc_6)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Conc_6)<-studentsname


for (f in 1:numstudents){
  matrix.Conc_6[f,1]<-mean(Conc_6[f,])
  matrix.Conc_6[f,2]<-sd(Conc_6[f,])
  matrix.Conc_6[f,3]<-min(Conc_6[f,])
  matrix.Conc_6[f,4]<-max(Conc_6[f,])
}

#Dose
matrix.Dose_6<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Dose_6)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Dose_6)<-studentsname

for (f in 1:numstudents){
  matrix.Dose_6[f,1]<-mean(Dose_6[f,])
  matrix.Dose_6[f,2]<-sd(Dose_6[f,])
  matrix.Dose_6[f,3]<-min(Dose_6[f,])
  matrix.Dose_6[f,4]<-max(Dose_6[f,])
}

#Risk 
matrix.Risk_6<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Risk_6)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Risk_6)<-studentsname

for (f in 1:numstudents){
  matrix.Risk_6[f,1]<-mean(Risk_6[f,])
  matrix.Risk_6[f,2]<-sd(Risk_6[f,])
  matrix.Risk_6[f,3]<-min(Risk_6[f,])
  matrix.Risk_6[f,4]<-max(Risk_6[f,])
}
#Check the data 
View(matrix.Conc_6)
View(matrix.Dose_6)
View(matrix.Risk_6)

#Pull out the data 
library(openxlsx)
write.csv(matrix.Conc_6, file="Conc_6.csv")
write.csv(matrix.Dose_6, file="Dose_6.csv")
write.csv(matrix.Risk_6, file="Risk_6.csv")

##6.4 8ACH data------------------------------------------- 

#Conc
matrix.Conc_8<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Conc_8)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Conc_8)<-studentsname


for (f in 1:numstudents){
  matrix.Conc_8[f,1]<-mean(Conc_8[f,])
  matrix.Conc_8[f,2]<-sd(Conc_8[f,])
  matrix.Conc_8[f,3]<-min(Conc_8[f,])
  matrix.Conc_8[f,4]<-max(Conc_8[f,])
}

#Dose
matrix.Dose_8<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Dose_8)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Dose_8)<-studentsname

for (f in 1:numstudents){
  matrix.Dose_8[f,1]<-mean(Dose_8[f,])
  matrix.Dose_8[f,2]<-sd(Dose_8[f,])
  matrix.Dose_8[f,3]<-min(Dose_8[f,])
  matrix.Dose_8[f,4]<-max(Dose_8[f,])
}

#Risk 
matrix.Risk_8<-matrix(nrow=numstudents, ncol=4)
colnames(matrix.Risk_8)<-c('mean', 'sd', 'min', 'max')
rownames(matrix.Risk_8)<-studentsname

for (f in 1:numstudents){
  matrix.Risk_8[f,1]<-mean(Risk_8[f,])
  matrix.Risk_8[f,2]<-sd(Risk_8[f,])
  matrix.Risk_8[f,3]<-min(Risk_8[f,])
  matrix.Risk_8[f,4]<-max(Risk_8[f,])
}
#Check the data 
View(matrix.Conc_8)
View(matrix.Dose_8)
View(matrix.Risk_8)

#Pull out the data 
library(openxlsx)
write.csv(matrix.Conc_8, file="Conc_8.csv")
write.csv(matrix.Dose_8, file="Dose_8.csv")
write.csv(matrix.Risk_8, file="Risk_8.csv")


#7. Sensitivity Analysis (Need to revise for this scenario)---------------------------------------------------

spear.Ecol<-data.frame(T.handarea, Frac.HS, Frac.HF, Reduc.intv,TE.all, TE.face,
                       Conc.i.face, Conc.i.hand, Conc.i.surface, Risk[2,])  

spear.anal<-cor(spear.Ecol,method="spearman")

View(spear.anal)

library(openxlsx)
write.csv (spear.anal, file="Sensitivity.ecol.csv")

