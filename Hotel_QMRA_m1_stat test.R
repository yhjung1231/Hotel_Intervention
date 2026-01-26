source("Hotel_Parameters_ver2.R")

# Touch scenario function: loop by the sequence and applying intervention factors

## sequence: character vector, e.g. c("Elevator"), c("Elevator", "Frontdesk")
## scenario: "baseline", "I1", "I2", "I3", "I4", "I5"
## Baseline - no intervention 
## I1: Using Antimicrobial film on the elevator button (Engineering control)
## I2: Targeted Hygiene on the elevator button (Administrative 1)
## I3: Using hand gloves - susceptible only (PPE)
## I4, 6, 8: Using hand sanitizer after touching elevator button - susceptible guest (Administrative 2)
## I5, 7, 9: Using hand sanitizer before touching elevator button - infected guest (Engineering 2- source control)

### I4,5 - P_HS:0.12/ I6,7 - P_HS:0.2/ I8,9- P_HS: 0.48


# Step 1: Main function for touch sequence 
run_touch_sequence <- function(sequence, scenario, iter, intervention = NULL) {
  
  numevents <- length(sequence) + 1  # +1 for face touch at the end
  eventsname <- c(paste0("suscept touch seq ", seq_along(sequence)), "Hand to face touch")
  
  Conc.h <- matrix(NA, nrow=numevents, ncol=iter, dimnames = list(eventsname, NULL))
  Conc.s <- matrix(NA, nrow=numevents, ncol=iter, dimnames = list(eventsname, NULL))
  Dose <- matrix(NA, nrow=numevents, ncol=iter, dimnames = list(eventsname, NULL))
  Risk <- matrix(NA, nrow=numevents, ncol=iter, dimnames = list(eventsname, NULL))
  
  #-------------------------------------------------------
  # Scene 1: Event 0 - Infected guest touches first surface
  #-------------------------------------------------------
  surf1 <- "Elevator"
  param1 <- get_surface_params(surf1, iter)
  
  # Transfer efficiencies (adjusted by intervention)
  TE.sh1 <- if (scenario == "I1") TE.h_fil 
  else param1$TE.sh
  
  TE.hs1 <- if (scenario == "I1") TE.fil_h 
  else param1$TE.hs
  
  # Surface inactivation rate
  k.surf <- if (scenario == "I1") k.surf.cu else param1$k.surf
  
  # Intervention 2: Targeted Hygiene (Administrative)
  LR.factor <- if (scenario == "I2") LR.S else 1
  
  # Initial surface concentration by infected person
  Conc.surf.i <- TE.hs1 * param1$Frac.hs * (T.handarea / param1$T.surfarea) * Conc.h.inf / LR.factor
  
  #Intervention 5, 7, 9: Source control (Engineering)
  P_HS <- if (scenario %in% c("I4", "I5")) 0.12 else
    if (scenario %in% c("I6", "I7")) 0.20 else
      if (scenario %in% c("I8", "I9")) 0.48 else 0
  
  
  if (scenario %in% c("I5", "I7", "I9")) {
    compliance <-runif(iter) < P_HS
    LR_applied <-ifelse (compliance, LR.HS, 1)
    Conc.surf.i <- Conc.surf.i / LR_applied
  }
  
  #----------------------------------------------------------------------------
  # Scene 2: Susceptible person touches surfaces (Loop over each touch sequence)
  #----------------------------------------------------------------------------
  # Initialize hand concentration before first surface contact
  # This is only used for the FIRST contact (sequence[1])
  # It will be updated with each loop iteration using Conc.h[i, ]
  
  Conc.h.prev <- Conc.h.sus
  
  for (i in seq_along(sequence)) {
    surf <- sequence[i]
    param <- get_surface_params(surf, iter)
    k.surf <- if (scenario == "I1" & surf == "Elevator") k.surf.cu else param$k.surf
    
    # Transfer efficiencies
    TE.sh <- if (scenario == "I3") param$TE.sh * hand_glove else if (scenario == "I1" & surf == "Elevator") TE.h_fil else param$TE.sh
    TE.hs <- if (scenario == "I3") param$TE.hs * hand_glove else if (scenario == "I1" & surf == "Elevator") TE.fil_h else param$TE.hs
    
    #Saving parameters
    debug_log_list[[length(debug_log_list) + 1]] <<- data.frame(
      Scenario = scenario,
      Surface = surf,
      TE.sh = mean(TE.sh),
      TE.hs = mean(TE.hs),
      k.surf = mean(k.surf),
      Conc.surf.i = mean(Conc.surf.i)
    )
    
    
    # Hand concentration after touching surface
    Conc.h[i, ] <- Conc.h.prev * exp(-k.hand * Time.m1) - {
      param$Frac.hs * (TE.hs * Conc.h.prev * exp(-k.hand * Time.m1) -
                         TE.sh * if (i == 1) Conc.surf.i * exp(-k.surf * Time.m1) else 0)
    }
    
    # Intervention 4, 6, 8: Hand hygiene after first elevator contact (Administrative)
    
    if (scenario %in% c("I4", "I6", "I8") && i == 1 && surf == "Elevator") {
      compliance <-runif(iter) < P_HS
      LR_applied <-ifelse(compliance, LR.HS, 1)
      Conc.h[i, ] <- Conc.h[i, ] / LR_applied
    }
    
    # Surface concentration update
    Conc.s[i, ] <- if (i == 1) Conc.surf.i * exp(-k.surf * Time.m1) else 0
    Conc.s[i, ] <- Conc.s[i, ] - {
      param$Frac.hs * (T.handarea / param$T.surfarea) *
        (TE.sh * if (i == 1) Conc.surf.i * exp(-k.surf * Time.m1) else 0 -
           Conc.h.prev * exp(-k.hand * Time.m1))
    }
    
    Dose[i, ] <- 0
    Risk[i, ] <- 0
    
    # Update previous hand concentration for next contact event
    Conc.h.prev <- Conc.h[i, ]
  }
  
  #------------------------------------
  # Scene 3: Final hand-to-face contact
  #------------------------------------
  Conc.h[numevents, ] <- (1 - TE.hf * Frac.hf) * (Conc.h.prev * exp(-k.hand * Time.m1))
  Conc.s[numevents, ] <- Conc.s[numevents - 1, ] * exp(-k.surf * Time.m1)
  Dose[numevents, ] <- TE.hf * Frac.hf * T.handarea * Conc.h.prev * exp(-k.hand * Time.m1)
  Risk[numevents, ] <- safe_risk_calc(Dose[numevents, ])
  
  return(list(Conc.h = Conc.h, Conc.s = Conc.s, Dose = Dose, Risk = Risk))
}

#=============================================================================================
# Step 2: Model fucntion to pack and run all the scenario in different sequence (sub-scenarios)
#==============================================================================================
run_model1_all <- function(iter) {
  
  # scenarios: baseline + I1~9
  scenarios <- c("Baseline", "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9")
  
  # sequence list 
  sequences <- list(
    E = c("Elevator"),
    E_FD = c("Elevator", "Frontdesk"),
    E_TT = c("Elevator", "Table")
  )
  
  results <- list()
  
  for (sc in scenarios) {
    for (seq_name in names(sequences)) {
      seq_val <- sequences[[seq_name]]
      out <- run_touch_sequence(seq_val, sc, iter)
      key <- paste0("Model1_", sc, "_", seq_name)
      results[[key]] <- out
    }
  }
  
  return(results)
}


#=====================================================
# Step 3: Prepare fuctions for plotting infection risk 
#=====================================================

library(ggplot2)
library(dplyr)
library(tidyr)

prepare_risk_plot_data <- function(results) {
  df_list <- list()
  
  for (name in names(results)) {
    res <- results[[name]]
    scenario_seq <- strsplit(name, "_")[[1]]  # e.g., Model1_I2_E_FD
    scenario <- scenario_seq[2]
    sequence <- paste(scenario_seq[3:length(scenario_seq)], collapse = "_")
    
    # Extract Risk values for the final hand-to-face event
    risk <- res$Risk["Hand to face touch", ]
    df <- data.frame(
      Scenario = scenario,
      Sequence = sequence,
      Risk = risk
    )
    df_list[[name]] <- df
  }
  
  # Combine all dataframes into one dataframe
  full_df <- do.call(rbind, df_list)
  return(full_df)
}


#========================================
#Step 4: Run model 1 
#========================================

# Run model
results <- run_model1_all(iter = 10000)


# Data prep for Stat test

risk_df <- prepare_risk_plot_data(results)

#========================================
#Step 5: Statistical Test  
#========================================
 #Difference between touch sequence scenarios 

library(dplyr)

kruskal.test(Risk ~ Sequence,
             data = risk_df %>% filter(Scenario == "Baseline"))

scenarios_to_check <- c("I1","I2","I3")

sequence_tests <- lapply(scenarios_to_check, function(sc) {
  out <- kruskal.test(
    Risk ~ Sequence,
    data = risk_df %>% dplyr::filter(Scenario == sc)
  )
  
  data.frame(
    Scenario = sc,          # == 말고 = !
    p_value  = out$p.value  # == 말고 = !
  )
})

sequence_tests_df <- do.call(rbind, sequence_tests)
sequence_tests_df


#Difference between hand hygien compliance rate scenario 

sequences_to_check<-c("E", "E_FD", "E_TT")


I468_tests <- lapply(sequences_to_check, function(seq) {
  out <- kruskal.test(
    Risk ~ Scenario,
    data = risk_df %>%
      filter(Sequence == seq, Scenario %in% c("I4","I6","I8"))
  )
  data.frame(
    Sequence = seq,
    p_value  = out$p.value
  )
})

I468_tests_df <- do.call(rbind, I468_tests)
I468_tests_df

library(FSA)

dunnTest(
  Risk ~ Scenario,
  data = risk_df %>%
    filter(Sequence == "E", Scenario %in% c("I4","I6","I8")),
  method = "holm"
)

dunnTest(
  Risk ~ Scenario,
  data = risk_df %>%
    filter(Sequence == "E_FD", Scenario %in% c("I4","I6","I8")),
  method = "holm"
)

dunnTest(
  Risk ~ Scenario,
  data = risk_df %>%
    filter(Sequence == "E_TT", Scenario %in% c("I4","I6","I8")),
  method = "holm"
)



I579_tests <- lapply(sequences_to_check, function(seq) {
  out <- kruskal.test(
    Risk ~ Scenario,
    data = risk_df %>%
      filter(Sequence == seq, Scenario %in% c("I5","I7","I9"))
  )
  data.frame(
    Sequence = seq,
    p_value  = out$p.value
  )
})

I579_tests_df <- do.call(rbind, I579_tests)
I579_tests_df




