#File 1 Parameters and distributions 
library(readxl)
library (truncdist)
library(triangle)
library(EnvStats)

iter <- 10000
set.seed(456)

#Norovirus concentration on hand of infected individual (genome copies per hand)
Conc.h.inf_gc<-10^(runif(iter, min=2.4, max=7.9)) # because log normal distribution is usually assumed for microorganisms.

#Recovery rate from sampling by sampling methods
RE.rinse<-runif(iter, 0.2, 0.84)
RE.swab<-runif(iter, 0.28,0.61)

#Ratio gc to PFU 
#gc_PFU<-10^(runif(iter,3.66, 4.46))
gc_PFU<-rtri(iter, 5.8, 47, 28)

#Total hand area (cm^2)
T.handarea<-runif(iter, min=445, max=535) 

#Norovirus concentration on hand of infected individual  (PFU / cm2)
Conc.h.inf <-Conc.h.inf_gc/(T.handarea*gc_PFU*RE.rinse)

#Norovirus concentration on the seeded surfaces (PFU/cm2)
Conc.seed<-1.87E+07

#Fraction of hand surface(hand-face contact) (unitless)
Frac.hf <-runif(iter, min=0.008, max=0.012) 

#Transfer efficiency (Hand to face)
TE.hf<-rtrunc(iter,"norm", mean=0.3390, sd=0.1318, a=0, b=1)

#Initial concentration on hand/face of susceptible person
Conc.h.sus<-0
Conc.f.sus<-0

#Inactivation rate on hands (k, min-1) 
x=2
Time.inact.h<- runif(iter, min=20, max=60)
k.hand<-log(x)/Time.inact.h

#Time between event 
Time.m1<-3

#Safe Risk calculation (Dose-response)

safe_risk_calc<-function(dose, alpha=0.722, beta=1106) {
  ratio <-dose/beta
  #Use linear approximation for very small ratio values: 1-exp(-x) ≈ x when x<<1
  ifelse (ratio < 1e-10,
          alpha*ratio, # linear approximation
          alpha*(1-exp(-ratio))) # accurate calculation
}


###Interventions---------------------------------------

#Engineering (Using Antimicrobial film)
k.surf.cu<-runif(iter, 0.2, 0.4)
TE.fil_h <-rtrunc(iter, "norm", mean=0.34, sd=0.12, a=0, b=1)
TE.h_fil <-rtrunc(iter, "norm", mean=0.13, sd=0.12, a=0, b=1)


#Administrative (Targeted Hygiene)
LR.S <- 10^(1.58)

# Administrative2: Hand Sanitizer
LR.HS <- 10^(rtrunc(iter, "norm", mean=1.06, sd=0.54, a=0.15, b=1.89))
P_HW<-0.3 # Hand hygiene compliance 

#PPE (gloves)
hand_glove<- runif(iter, 0.53,0.68)

TE.gf<-TE.hf*hand_glove

t_gl=3*24*60

k_gl=-log(10^(-3.5))/t_gl


#_________________________________________________________


## Surface Specific Parameters (Updated 052325)
if (fomite == "Elevator"){
  
  #Surface area (cm^2)
  T.surfarea<-10.14 
  
  #Fraction of hand surface (hand-surface) (unitless)
  Frac.hs <-runif(iter, min=0.008, max=0.012) 
  
  #Recovered PhiX 174 Concentration after 4 hour (PFU/cm2)
  Conc.recover<-2.47
  
  #Transfer efficiency (surface to hand, hand to surface) (unitless)
  TE.sh<-rtrunc(iter,"norm", mean=0.37, sd=0.14, a=0, b=1)
  TE.hs<-rtrunc(iter,"norm", mean=0.16, sd=0.16, a=0, b=1)
  
  TE.sg<-TE.sh*hand_glove
  TE.gs<-TE.hs*hand_glove
  
  #Time of inactivation on surface i by 1/x amount of initial concentration (min)
  x=2
  Time.inact<- 2*24*60
 
  #iactivation rate (k, min-1)
  k.surf<-log(x)/Time.inact
  
} else if (fomite == "Frontdesk"){
  
  #Surface area (cm^2)
  T.surfarea<-871.2 
  
  #Fraction of hand surface (hand-surface) (unitless)
  Frac.hs <-runif(iter, min=0.013, max=0.25) 
  
  #Recovered PhiX 174 Concentration after 4 hour (PFU/cm2)
  Conc.recover<-1.15E-01
  
  #Transfer efficiency (surface to hand, hand to surface) (unitless)
  TE.sh<-rtrunc(iter,"norm", mean=0.37, sd=0.14, a=0, b=1)
  TE.hs<-rtrunc(iter,"norm", mean=0.16, sd=0.16, a=0, b=1)
 
  TE.sg<-TE.sh*hand_glove
  TE.gs<-TE.hs*hand_glove
  
  
  #Time of inactivation on surface i by 1/x amount of initial concentration (min)
  x=2
  Time.inact<- 2*24*60
  #iactivation rate (k, min-1)
  k.surf<-log(x)/Time.inact
  
  
}else {#fomite == "Table"
  
  #Surface area (cm^2)
  T.surfarea<-runif(iter, 4879, 12053)
  
  #Fraction of hand surface (hand-surface) (unitless)
  Frac.hs <-runif(iter, min=0.013, max=0.25) 
  
  #Recovered PhiX 174 Concentration after 4 hour (PFU/cm2)
  Conc.recover<-runif(iter, 1.04E-04, 1.08E-03)
  
  #Transfer efficiency (surface to hand, hand to surface) (unitless)
  TE.sh<-rtrunc(iter,"norm", mean=0.3, sd=0.18, a=0, b=1)
  TE.hs<-rtrunc(iter,"norm", mean=0.22, sd=0.17, a=0, b=1)
  
  TE.sg<-TE.sh*hand_glove
  TE.gs<-TE.hs*hand_glove
  
  #Time of inactivation on surface i by 1/x amount of initial concentration (min)
  x=10
  Time.inact<- 7.2*24*60
  #iactivation rate (k, min-1)
  k.surf<-log(x)/Time.inact
}





