# File 1: Parameters and Distributions ---------------------------------------

if (!require(pacman)) install.packages("pacman")
pacman::p_load(readxl, truncdist, triangle,EnvStats)

iter <- 10000
set.seed(456)

# Norovirus concentration on hand of infected individual (gc/hand)
Conc.h.inf_gc <- 10^(runif(iter, min=2.4, max=7.9))  # log-normal assumption

# Recovery efficiency
RE.rinse <- runif(iter, 0.2, 0.84)

# Ratio gc to PFU
gc_PFU <- rtri(iter, 5.8, 47, 28)

# Total hand area (cm^2)
T.handarea <- runif(iter, 445, 535)

# Norovirus concentration on hand (PFU/cm2)
Conc.h.inf <- Conc.h.inf_gc / (T.handarea * gc_PFU * RE.rinse)

# Seeded surface concentration (PFU/cm2)
Conc.samp.elv <- rtri(iter, 5.0E+05, 5.0E+06, 1.0E+06)
Conc.seed <- Conc.samp.elv/10.14 

# Fraction of hand touching face
Frac.hf <- runif(iter, 0.008, 0.012)

# Transfer efficiency: hand to face
TE.hf <- rtrunc(iter, "norm", mean=0.3390, sd=0.1318, a=0, b=1)

# Initial concentrations for susceptible person
Conc.h.sus <- 0
Conc.f.sus <- 0

# Inactivation rate on hand
x <- 2
Time.inact.h <- runif(iter, 20, 60)
k.hand <- log(x) / Time.inact.h

# Time between events (min)
Time.m1 <- 3

# Dose-response function
safe_risk_calc <- function(dose, alpha=0.722, beta=1106) {
  ratio <- dose / beta
  ifelse(ratio < 1e-10,
         alpha * ratio,
         alpha * (1 - exp(-ratio)))
}

# Interventions --------------------------------------------------------------

# Engineering: Antimicrobial surface (copper)
k.surf.cu <- runif(iter, 0.2, 0.4)
TE.fil_h <-rtrunc(iter, "norm", mean=0.34, sd=0.12, a=0, b=1)
TE.h_fil <-rtrunc(iter, "norm", mean=0.13, sd=0.12, a=0, b=1)

# Administrative1: Targeted hygiene (log reduction)
LR.S <- 10^(rtri(iter, 5.93E-03, 2.57, 1.58))  

# Administrative2: Hand Sanitizer
LR.HS <- 10^(rtrunc(iter, "norm", mean=1.06, sd=0.54, a=0.15, b=1.89))
## Hand hygiene compliance (P_HW) would be added in the model  

# PPE: Gloves (TE.surf of glove are existed under the get_surface_params function)
hand_glove <- runif(iter, 0.53, 0.68)  # effectiveness of gloves
TE.gf <- TE.hf * hand_glove           # hand-to-face with gloves

# Inactivation rate on glove
t_gl <- runif(iter,2,4)*24 * 60 #2-4 days
k_gl <- -log(10^(-3.5)) / t_gl

# Surface-specific parameters function ---------------------------------------

get_surface_params <- function(fomite, iter) {
  Sample.surfarea <- NA 
  if (fomite == "Elevator") {
    T.surfarea <- 10.14
    Sample.surfarea<-T.surfarea
    Frac.hs <- runif(iter, 0.008, 0.012)
    Conc.recover <- runif(iter, 0, 25/Sample.surfarea)
    TE.sh <- rtrunc(iter, "norm", mean=0.37, sd=0.14, a=0, b=1)
    TE.hs <- rtrunc(iter, "norm", mean=0.16, sd=0.16, a=0, b=1)
    x <- 2
    Time.inact <- runif(iter, 38, 58) * 60
    
  } else if (fomite == "Frontdesk") {
    T.surfarea <- 871.2
    Sample.surfarea<-T.surfarea/4
    Frac.hs <- runif(iter, 0.013, 0.25)
    Conc.recover <-runif(iter, 0, 25/Sample.surfarea)
    TE.sh <- rtrunc(iter, "norm", mean=0.37, sd=0.14, a=0, b=1)
    TE.hs <- rtrunc(iter, "norm", mean=0.16, sd=0.16, a=0, b=1)
    x <- 2
    Time.inact <- runif(iter, 38, 58) * 60
    
  } else if (fomite == "Table") {
    T.surfarea <- runif(iter, 4879, 14457)
    Sample.surfarea <-T.surfarea/4
    Frac.hs <- runif(iter, 0.013, 0.25)
    Conc.recover <- rtri(iter, 0, 4.31E-03, 5.11E-04)
    TE.sh <- rtrunc(iter, "norm", mean=0.3, sd=0.18, a=0, b=1)
    TE.hs <- rtrunc(iter, "norm", mean=0.22, sd=0.17, a=0, b=1)
    x <- 10
    Time.inact <- runif(iter, 6.2, 8.2) * 24 * 60
    
  } else {
    stop("Invalid fomite type: should be 'Elevator', 'Frontdesk', or 'Table'")
  }
  
  k.surf <- log(x) / Time.inact
  
  return(list(
    T.surfarea = T.surfarea,
    Sample.surfarea = Sample.surfarea,
    Frac.hs = Frac.hs,
    Conc.recover = Conc.recover,
    TE.sh = TE.sh,
    TE.hs = TE.hs,
    k.surf = k.surf
  ))
}

# Example of how to use the parameter function -------------------------------

#fomite <- "Elevator"  # or "Frontdesk", "Table"
#params <- get_surface_params(fomite, iter)

# Derive glove-based transfer efficiencies (outside the function)
TE.sg <- params$TE.sh * hand_glove
TE.gs <- params$TE.hs * hand_glove

# Now you can use:
# - params$T.surfarea
# - params$Frac.hs
# - params$Conc.recover
# - params$k.surf
# - TE.sg, TE.gs

#summary of starting concentration info 
params_E<- get_surface_params ("Elevator", iter)
params_FD<-get_surface_params ("Frontdesk", iter)
params_TT<-get_surface_params ("Table", iter)

#M1
summarise(params_E$Conc.surf.i)
