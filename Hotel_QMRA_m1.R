fomite <-"Elevator"
source("Hotel_Parameters.R")

###Model 1: Susceptible touch surfaces right after an infected guest passes by the hotel lobby###

## Baseline scenario=========================================================================

##M1.S1: "Elevator" touch____________________________________________   

###S1.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.1<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.1)<-eventsname

Conc.s.1<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.1)<-eventsname

Dose.1<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.1)<-eventsname

Risk.1<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.1)<-eventsname

#Event 0. Infected person passed by hotel lobby

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf

#Event 1. Susceptible guest contact elevator button right after infected guest passed by the lobby area
Conc.h.1[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.surf.i*exp(-k.surf*Time.m1))}
Conc.s.1[1,]<-Conc.surf.i*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.surf.i*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.1[1,]<-0
Risk.1[1,]<-0

#Event 2. Hand to face contact
Conc.h.1[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.1[1,]*exp(-k.hand*Time.m1))
Conc.s.1[2,]<-Conc.s.1[1,]*exp(-k.surf*Time.m1)
Dose.1[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.1[1,]*exp(-k.hand*Time.m1)
Risk.1[2,]<-safe_risk_calc(Dose.1[2,])


##M1.S2: "Front desk" touch/ ##M1.S4: "Table top" touch__________________________________

###s1.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.2<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.2)<-eventsname

Conc.s.2<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.2)<-eventsname

Dose.2<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.2)<-eventsname

Risk.2<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.2)<-eventsname

#Event 0. Infected person passed by hotel lobby

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf


#Event 1. Susceptible guest contact Front desk right after infected guest passed by the lobby area
fomite <-"Frontdesk"
source("Hotel_Parameters.R")

Conc.h.2[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.2[1,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.2[1,]<-0
Risk.2[1,]<-0

#Event 2. Hand to face contact
Conc.h.2[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.2[1,]*exp(-k.hand*Time.m1))
Conc.s.2[2,]<-Conc.s.2[1,]*exp(-k.surf*Time.m1)
Dose.2[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.2[1,]*exp(-k.hand*Time.m1)
Risk.2[2,]<-safe_risk_calc(Dose.2[2,])

##M1.S3: "Elevator - Front desk" touch__________________________________
###s3.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.3<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.3)<-eventsname

Conc.s.3<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.3)<-eventsname

Dose.3<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.3)<-eventsname

Risk.3<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.3)<-eventsname

#Event 0. Infected person passed by hotel lobby
fomite <-"Elevator"
source("Hotel_Parameters.R")

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf


#Event 1. Susceptible guest contact seq 1 (Elevator)
Conc.h.3[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.surf.i*exp(-k.surf*Time.m1))}
Conc.s.3[1,]<-Conc.surf.i*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.surf.i*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.3[1,]<-0
Risk.3[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Front desk)
fomite <-"Frontdesk"
source("Hotel_Parameters.R")

Conc.h.3[2,]<-Conc.h.3[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.3[1,]*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.3[2,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.3[1,]*exp(-k.hand*Time.m1))}
Dose.3[2,]<-0
Risk.3[2,]<-0


#Event 3. Hand to face contact
Conc.h.3[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.3[2,]*exp(-k.hand*Time.m1))
Conc.s.3[3,]<-Conc.s.3[2,]*exp(-k.surf*Time.m1)
Dose.3[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.3[2,]*exp(-k.hand*Time.m1)
Risk.3[3,]<-safe_risk_calc(Dose.3[3,])

##M1.S5: "Elevator - Table top" touch___________________________________

###s5.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.5<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.5)<-eventsname

Conc.s.5<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.5)<-eventsname

Dose.5<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.5)<-eventsname

Risk.5<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.5)<-eventsname

#Event 0. Infected person passed by hotel lobby
fomite <-"Elevator"
source("Hotel_Parameters.R")

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf


#Event 1. Susceptible guest contact seq 1 (Elevator)
Conc.h.5[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.surf.i*exp(-k.surf*Time.m1))}
Conc.s.5[1,]<-Conc.surf.i*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.surf.i*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.5[1,]<-0
Risk.5[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Table))
fomite <-"Table"
source("Hotel_Parameters.R")

Conc.h.5[2,]<-Conc.h.5[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.5[1,]*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.5[2,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.5[1,]*exp(-k.hand*Time.m1))}
Dose.5[2,]<-0
Risk.5[2,]<-0


#Event 3. Hand to face contact
Conc.h.5[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.5[2,]*exp(-k.hand*Time.m1))
Conc.s.5[3,]<-Conc.s.5[2,]*exp(-k.surf*Time.m1)
Dose.5[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.5[2,]*exp(-k.hand*Time.m1)
Risk.5[3,]<-safe_risk_calc(Dose.5[3,])

## Intervention1: Engineering scenario===============================================

#change k.surf parameter on elevator to k.surf.cu
#Change TE.hs, TE.sh --> TE.h_fil, TE.fil_h

##M1.S1.I1: "Elevator" touch____________________________________________   

###S1.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.1i<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.1i)<-eventsname

Conc.s.1i<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.1i)<-eventsname

Dose.1i<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.1i)<-eventsname

Risk.1i<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.1i)<-eventsname

#Event 0. Infected person passed by hotel lobby

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.h_fil*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf

#Event 1. Susceptible guest contact elevator button right after infected guest passed by the lobby area
Conc.h.1i[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.h_fil*Conc.h.sus*exp(-k.hand*Time.m1)-TE.fil_h*Conc.surf.i*exp(-k.surf.cu*Time.m1))}
Conc.s.1i[1,]<-Conc.surf.i*exp(-k.surf.cu*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.fil_h*Conc.surf.i*exp(-k.surf.cu*Time.m1)-TE.h_fil*Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.1i[1,]<-0
Risk.1i[1,]<-0

#Event 2. Hand to face contact
Conc.h.1i[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.1i[1,]*exp(-k.hand*Time.m1))
Conc.s.1i[2,]<-Conc.s.1i[1,]*exp(-k.surf.cu*Time.m1)
Dose.1i[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.1i[1,]*exp(-k.hand*Time.m1)
Risk.1i[2,]<-safe_risk_calc(Dose.1i[2,])



##M1.S3.I1: "Elevator - Front desk" touch__________________________________
###s3.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.3i<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.3i)<-eventsname

Conc.s.3i<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.3i)<-eventsname

Dose.3i<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.3i)<-eventsname

Risk.3i<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.3i)<-eventsname

#Event 0. Infected person passed by hotel lobby
fomite <-"Elevator"
source("Hotel_Parameters.R")

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.h_fil*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf

# k.surf.cu was applied only for elevator button

#Event 1. Susceptible guest contact seq 1 (Elevator)
Conc.h.3i[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.h_fil*Conc.h.sus*exp(-k.hand*Time.m1)-TE.fil_h*Conc.surf.i*exp(-k.surf.cu*Time.m1))}
Conc.s.3i[1,]<-Conc.surf.i*exp(-k.surf.cu*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.fil_h*Conc.surf.i*exp(-k.surf.cu*Time.m1)-TE.h_fil*Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.3i[1,]<-0
Risk.3i[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Front desk)
fomite <-"Frontdesk"
source("Hotel_Parameters.R")

Conc.h.3i[2,]<-Conc.h.3i[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.3i[1,]*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.3i[2,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.3i[1,]*exp(-k.hand*Time.m1))}
Dose.3i[2,]<-0
Risk.3i[2,]<-0


#Event 3. Hand to face contact
Conc.h.3i[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.3i[2,]*exp(-k.hand*Time.m1))
Conc.s.3i[3,]<-Conc.s.3i[2,]*exp(-k.surf*Time.m1)
Dose.3i[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.3i[2,]*exp(-k.hand*Time.m1)
Risk.3i[3,]<-safe_risk_calc(Dose.3i[3,])


##M1.S5.I1: "Elevator - Table top" touch___________________________________

###s5.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.5i<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.5i)<-eventsname

Conc.s.5i<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.5i)<-eventsname

Dose.5i<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.5i)<-eventsname

Risk.5i<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.5i)<-eventsname

#Event 0. Infected person passed by hotel lobby
fomite <-"Elevator"
source("Hotel_Parameters.R")

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.h_fil*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf

# k.surf.cu and TE.h_fil/TE.fil_h were applied only for elevator button

#Event 1. Susceptible guest contact seq 1 (Elevator)
Conc.h.5i[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.h_fil*Conc.h.sus*exp(-k.hand*Time.m1)-TE.fil_h*Conc.surf.i*exp(-k.surf.cu*Time.m1))}
Conc.s.5i[1,]<-Conc.surf.i*exp(-k.surf.cu*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.fil_h*Conc.surf.i*exp(-k.surf.cu*Time.m1)-TE.h_fil*Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.5i[1,]<-0
Risk.5i[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Table))
fomite <-"Table"
source("Hotel_Parameters.R")

Conc.h.5i[2,]<-Conc.h.5i[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.5i[1,]*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.5i[2,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.5i[1,]*exp(-k.hand*Time.m1))}
Dose.5i[2,]<-0
Risk.5i[2,]<-0


#Event 3. Hand to face contact
Conc.h.5i[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.5i[2,]*exp(-k.hand*Time.m1))
Conc.s.5i[3,]<-Conc.s.5i[2,]*exp(-k.surf*Time.m1)
Dose.5i[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.5i[2,]*exp(-k.hand*Time.m1)
Risk.5i[3,]<-safe_risk_calc(Dose.5i[3,])

## Intervention2: Administrative scenario============================(Updating 6/13/2025)
#Additional Targeted Hygiene was applied to Elevator button before susceptible hotel guest touch it.

##M1.S1.I2: "Elevator" touch____________________________________________   

###S1.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.1i2<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.1i2)<-eventsname

Conc.s.1i2<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.1i2)<-eventsname

Dose.1i2<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.1i2)<-eventsname

Risk.1i2<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.1i2)<-eventsname

#Event 0. Infected person passed by hotel lobby

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf/LR.S

#Event 1. Susceptible guest contact elevator button right after infected guest passed by the lobby area
Conc.h.11.1i2[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.surf.i*exp(-k.surf*Time.m1))}
Conc.s.1i2[1,]<-Conc.surf.i*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.surf.i*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.1i2[1,]<-0
Risk.1i2[1,]<-0

#Event 2. Hand to face contact
Conc.h.1i2[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.1i2[1,]*exp(-k.hand*Time.m1))
Conc.s.1i2[2,]<-Conc.s.1i2[1,]*exp(-k.surf*Time.m1)
Dose.1i2[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.1i2[1,]*exp(-k.hand*Time.m1)
Risk.1i2[2,]<-safe_risk_calc(Dose.1i2[2,])


##M1.S2.I2: "Front desk" touch/ ##M1.S4.I4: "Table top" touch__________________________________

###s1.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.2i2<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.2i2)<-eventsname

Conc.s.2i2<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.2i2)<-eventsname

Dose.2i2<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.2i2)<-eventsname

Risk.2i2<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.2i2)<-eventsname

#Event 0. Infected person passed by hotel lobby

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf/LR.S


#Event 1. Susceptible guest contact Front desk right after infected guest passed by the lobby area
fomite <-"Frontdesk"
source("Hotel_Parameters.R")

Conc.h.2i2[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.2i2[1,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.2i2[1,]<-0
Risk.2i2[1,]<-0

#Event 2. Hand to face contact
Conc.h.2i2[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.2i2[1,]*exp(-k.hand*Time.m1))
Conc.s.2i2[2,]<-Conc.s.2i2[1,]*exp(-k.surf*Time.m1)
Dose.2i2[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.2i2[1,]*exp(-k.hand*Time.m1)
Risk.2i2[2,]<-safe_risk_calc(Dose.2i2[2,])



##M1.S3.I2: "Elevator - Front desk" touch__________________________________
###s3.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.3i2<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.3i2)<-eventsname

Conc.s.3i2<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.3i2)<-eventsname

Dose.3i2<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.3i2)<-eventsname

Risk.3i2<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.3i2)<-eventsname

#Event 0. Infected person passed by hotel lobby
fomite <-"Elevator"
source("Hotel_Parameters.R")

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf/LR.S

#Event 1. Susceptible guest contact seq 1 (Elevator)
Conc.h.3i2[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.surf.i*exp(-k.surf*Time.m1))}
Conc.s.3i2[1,]<-Conc.surf.i*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.surf.i*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.3i2[1,]<-0
Risk.3i2[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Front desk)
fomite <-"Frontdesk"
source("Hotel_Parameters.R")

Conc.h.3i2[2,]<-Conc.h.3i2[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.3i2[1,]*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.3i2[2,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.3i2[1,]*exp(-k.hand*Time.m1))}
Dose.3i2[2,]<-0
Risk.3i2[2,]<-0


#Event 3. Hand to face contact
Conc.h.3i2[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.3i2[2,]*exp(-k.hand*Time.m1))
Conc.s.3i2[3,]<-Conc.s.3i2[2,]*exp(-k.surf*Time.m1)
Dose.3i2[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.3i2[2,]*exp(-k.hand*Time.m1)
Risk.3i2[3,]<-safe_risk_calc(Dose.3i2[3,])



##M1.S5.I2: "Elevator - Table top" touch___________________________________

###s5.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.5i2<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.5i2)<-eventsname

Conc.s.5i2<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.5i2)<-eventsname

Dose.5i2<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.5i2)<-eventsname

Risk.5i2<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.5i2)<-eventsname

#Event 0. Infected person passed by hotel lobby
fomite <-"Elevator"
source("Hotel_Parameters.R")

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf/LR.S


#Event 1. Susceptible guest contact seq 1 (Elevator)
Conc.h.5i2[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.surf.i*exp(-k.surf*Time.m1))}
Conc.s.5i2[1,]<-Conc.surf.i*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.surf.i*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.5i2[1,]<-0
Risk.5i2[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Table))
fomite <-"Table"
source("Hotel_Parameters.R")

Conc.h.5i2[2,]<-Conc.h.5i2[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.5i2[1,]*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.5i2[2,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.5i2[1,]*exp(-k.hand*Time.m1))}
Dose.5i2[2,]<-0
Risk.5i2[2,]<-0


#Event 3. Hand to face contact
Conc.h.5i2[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.5i2[2,]*exp(-k.hand*Time.m1))
Conc.s.5i2[3,]<-Conc.s.5i2[2,]*exp(-k.surf*Time.m1)
Dose.5i2[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.5i2[2,]*exp(-k.hand*Time.m1)
Risk.5i2[3,]<-safe_risk_calc(Dose.5i2[3,])

## Intervention3: PPE============================(Updating 6/13/2025)
#Additional Targeted Hygiene was applied to Elevator button before susceptible hotel guest touch it.

##M1.S1.I2: "Elevator" touch____________________________________________   

###S1.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.1i2<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.1i2)<-eventsname

Conc.s.1i2<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.1i2)<-eventsname

Dose.1i2<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.1i2)<-eventsname

Risk.1i2<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.1i2)<-eventsname

#Event 0. Infected person passed by hotel lobby

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf/LR.S

#Event 1. Susceptible guest contact elevator button right after infected guest passed by the lobby area
Conc.h.11.1i2[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.surf.i*exp(-k.surf*Time.m1))}
Conc.s.1i2[1,]<-Conc.surf.i*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.surf.i*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.1i2[1,]<-0
Risk.1i2[1,]<-0

#Event 2. Hand to face contact
Conc.h.1i2[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.1i2[1,]*exp(-k.hand*Time.m1))
Conc.s.1i2[2,]<-Conc.s.1i2[1,]*exp(-k.surf*Time.m1)
Dose.1i2[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.1i2[1,]*exp(-k.hand*Time.m1)
Risk.1i2[2,]<-safe_risk_calc(Dose.1i2[2,])


##M1.S2.I2: "Front desk" touch/ ##M1.S4.I4: "Table top" touch__________________________________

###s1.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","Hand to face touch")
numevents<-length(eventsname)

Conc.h.2i2<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.2i2)<-eventsname

Conc.s.2i2<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.2i2)<-eventsname

Dose.2i2<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.2i2)<-eventsname

Risk.2i2<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.2i2)<-eventsname

#Event 0. Infected person passed by hotel lobby

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf/LR.S


#Event 1. Susceptible guest contact Front desk right after infected guest passed by the lobby area
fomite <-"Frontdesk"
source("Hotel_Parameters.R")

Conc.h.2i2[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.2i2[1,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.2i2[1,]<-0
Risk.2i2[1,]<-0

#Event 2. Hand to face contact
Conc.h.2i2[2,]<-(1-TE.hf*Frac.hf)*(Conc.h.2i2[1,]*exp(-k.hand*Time.m1))
Conc.s.2i2[2,]<-Conc.s.2i2[1,]*exp(-k.surf*Time.m1)
Dose.2i2[2,]<-TE.hf*Frac.hf*T.handarea*Conc.h.2i2[1,]*exp(-k.hand*Time.m1)
Risk.2i2[2,]<-safe_risk_calc(Dose.2i2[2,])



##M1.S3.I3: "Elevator - Front desk" touch__________________________________
###s3.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.3i2<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.3i2)<-eventsname

Conc.s.3i2<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.3i2)<-eventsname

Dose.3i2<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.3i2)<-eventsname

Risk.3i2<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.3i2)<-eventsname

#Event 0. Infected person passed by hotel lobby
fomite <-"Elevator"
source("Hotel_Parameters.R")

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf/LR.S

#Event 1. Susceptible guest contact seq 1 (Elevator)
Conc.h.3i2[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.surf.i*exp(-k.surf*Time.m1))}
Conc.s.3i2[1,]<-Conc.surf.i*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.surf.i*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.3i2[1,]<-0
Risk.3i2[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Front desk)
fomite <-"Frontdesk"
source("Hotel_Parameters.R")

Conc.h.3i2[2,]<-Conc.h.3i2[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.3i2[1,]*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.3i2[2,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.3i2[1,]*exp(-k.hand*Time.m1))}
Dose.3i2[2,]<-0
Risk.3i2[2,]<-0


#Event 3. Hand to face contact
Conc.h.3i2[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.3i2[2,]*exp(-k.hand*Time.m1))
Conc.s.3i2[3,]<-Conc.s.3i2[2,]*exp(-k.surf*Time.m1)
Dose.3i2[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.3i2[2,]*exp(-k.hand*Time.m1)
Risk.3i2[3,]<-safe_risk_calc(Dose.3i2[3,])



##M1.S5.I5: "Elevator - Table top" touch___________________________________

###s5.0 Matrix -------------------------------------------------
eventsname<-c("suscept touch seq 1","suscept touch seq 2","Hand to face touch")
numevents<-length(eventsname)

Conc.h.5i2<-matrix(nrow=numevents,ncol=iter) # hand of susceptible guest
rownames(Conc.h.5i2)<-eventsname

Conc.s.5i2<-matrix(nrow=numevents, ncol=iter)
rownames(Conc.s.5i2)<-eventsname

Dose.5i2<-matrix(nrow=numevents, ncol=iter)
rownames(Dose.5i2)<-eventsname

Risk.5i2<-matrix(nrow=numevents, ncol=iter)
rownames(Risk.5i2)<-eventsname

#Event 0. Infected person passed by hotel lobby
fomite <-"Elevator"
source("Hotel_Parameters.R")

#Concentration of norovirus on surface which infected guest touched (Elevator& Back door)
Conc.surf.i<- TE.hs*Frac.hs*(T.handarea/T.surfarea)*Conc.h.inf/LR.S


#Event 1. Susceptible guest contact seq 1 (Elevator)
Conc.h.5i2[1,]<-Conc.h.sus*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.sus*exp(-k.hand*Time.m1)-TE.sh*Conc.surf.i*exp(-k.surf*Time.m1))}
Conc.s.5i2[1,]<-Conc.surf.i*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*Conc.surf.i*exp(-k.surf*Time.m1)-Conc.h.sus*exp(-k.hand*Time.m1))}
Dose.5i2[1,]<-0
Risk.5i2[1,]<-0

#Event 2. Susceptible guest contact seq 2 (Table))
fomite <-"Table"
source("Hotel_Parameters.R")

Conc.h.5i2[2,]<-Conc.h.5i2[1,]*exp(-k.hand*Time.m1)-{Frac.hs*(TE.hs*Conc.h.5i2[1,]*exp(-k.hand*Time.m1)-TE.sh*0*exp(-k.surf*Time.m1))}
Conc.s.5i2[2,]<-0*exp(-k.surf*Time.m1)-{Frac.hs*T.handarea/T.surfarea*(TE.sh*0*exp(-k.surf*Time.m1)-Conc.h.5i2[1,]*exp(-k.hand*Time.m1))}
Dose.5i2[2,]<-0
Risk.5i2[2,]<-0


#Event 3. Hand to face contact
Conc.h.5i2[3,]<-(1-TE.hf*Frac.hf)*(Conc.h.5i2[2,]*exp(-k.hand*Time.m1))
Conc.s.5i2[3,]<-Conc.s.5i2[2,]*exp(-k.surf*Time.m1)
Dose.5i2[3,]<-TE.hf*Frac.hf*T.handarea*Conc.h.5i2[2,]*exp(-k.hand*Time.m1)
Risk.5i2[3,]<-safe_risk_calc(Dose.5i2[3,])

##### Figure making and Result pulling ========================================================
##5. plotting--(have to update intervention result for plotting)----------------------------------------------------------
library(ggplot2)
library(ggpubr)


Risk.1.df<-as.data.frame(t(Risk.1))
Risk.3.df<-as.data.frame(t(Risk.3))
Risk.5.df<-as.data.frame(t(Risk.5))


type<-c(rep("E",iter),rep("E-FD",iter),rep("E-TT",iter))
value<-c(Risk.1.df$"Hand to face touch", Risk.3.df$"Hand to face touch", Risk.5.df$"Hand to face touch")  
        
data.m1<-data.frame(type,value)

windows()
ggplot(data.m1)+geom_violin(aes(x=type,y=value,fill=type, group=type),alpha=0.3,draw_quantiles = c(0.25,0.5,0.75))+
  scale_y_continuous(trans="log10", breaks=c( 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11),
                     label=c(-1:-11)) +
  scale_x_discrete(limits=c("E","E-FD","E-TT"))+
  ggtitle("Model 1: Susceptible guest touches fomite immediately after lobby use by infected guest")+
  xlab("Susceptible Guest Touch Sequences")+
  ylab(expression("Infection risk of norovirus (log" [10] * ")"))
  

ggsave("Model1.tiff", dpi=600, dev='tiff', height=6, width=8, units="in")



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

