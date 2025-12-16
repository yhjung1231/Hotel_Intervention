###======================================================================================###
#Model 2: Susceptible touch surfaces 4hour after an infected guest passes by the hotel lobby
###======================================================================================###

source("Hotel_Parameters_ver2.R")

# Touch scenario function: loop by the sequence and applying intervention factors

## sequence: character vector, e.g. c("Elevator"), c("Elevator", "Frontdesk")
## scenario: "Baseline", "I1", "I2", "I3", "I4", "I5"
## Baseline - no intervention 
## I1: Using Antimicrobial film on the elevator button (Engineering control)
## I2: Targeted Hygiene on the elevator button (Administrative 1)
## I3: Using hand gloves - susceptible only (PPE)
## I4: I1+I2
## I5: I1+I2+I3
## I6: I2+I3


debug_log_list<-list()

# Step 1: Main function for touch sequence 
run_touch_sequence_2 <- function(sequence, scenario, iter, intervention = NULL) {
  
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
  TE.sh1 <- if (scenario %in% c("I1", "I1+I2", "I1+I2+I3") && surf1 == "Elevator") TE.h_fil 
  else param1$TE.sh
  
  TE.hs1 <- if (scenario %in% c("I1", "I1+I2", "I1+I2+I3") && surf1 == "Elevator") TE.fil_h 
  else param1$TE.hs
  
  # Surface inactivation rate
  k.surf <- if (scenario %in% c("I1", "I1+I2", "I1+I2+I3") && surf1 == "Elevator") k.surf.cu else param1$k.surf
  
  # Intervention 2: Targeted Hygiene (Administrative)
  LR.factor <- if (scenario %in% c("I2", "I1+I2", "I1+I2+I3", "I2+I3") && surf1 == "Elevator") LR.S else 1
  
  # Initial surface concentration by infected person
  Conc.surf.i <- TE.hs1 * param1$Frac.hs * (T.handarea / param1$T.surfarea) * Conc.h.inf / LR.factor

  
  #----------------------------------------------------------------------------
  # Scene 2: Susceptible person touches surfaces (Loop over each touch sequence)
  #----------------------------------------------------------------------------
  # Initialize hand concentration before first surface contact
  # This is only used for the FIRST contact (sequence[1])
  # It will be updated with each loop iteration using Conc.h[i, ]
  
  Conc.h.prev <- Conc.h.sus
  Time.m2<-4*60 #(min, 4 hours)
  
  
  
  for (i in seq_along(sequence)) {
    surf <- sequence[i]
    param <- get_surface_params(surf, iter)
    
    
    k.surf <- if (scenario %in% c("I1", "I1+I2", "I1+I2+I3") && surf == "Elevator") k.surf.cu else param$k.surf
    Conc.recover <- param$Conc.recover
    
    
    #Concentration on each fomite after 4 hours 
    
    Conc.4hour<-if (scenario %in% c("I1", "I1+I2", "I1+I2+I3")&& surf == "Elevator"){(Conc.recover/Conc.seed)*Conc.surf.i*exp(-k.surf*Time.m2)}
    else {(Conc.recover/Conc.seed)*Conc.surf.i}
    
    
    # Transfer efficiencies
    TE.sh <- if (scenario %in% c("I3", "I2+I3")) param$TE.sh * hand_glove else if (scenario %in% c("I1", "I1+I2") & surf == "Elevator") TE.h_fil else if (scenario =="I1+I2+I3" & surf == "Elevator") TE.h_fil* hand_glove else param$TE.sh
    TE.hs <- if (scenario %in% c("I3", "I2+I3")) param$TE.hs * hand_glove else if (scenario %in% c("I1", "I1+I2") & surf == "Elevator") TE.fil_h else if (scenario =="I1+I2+I3" & surf == "Elevator") TE.h_fil* hand_glove else param$TE.hs
    
    #Debuging log / parameter checking
    debug_log_list[[length(debug_log_list) + 1]] <<- data.frame(
      Scenario = scenario,
      Surface = surf,
      TE.sh = mean(TE.sh),
      TE.hs = mean(TE.hs),
      k.surf = mean(k.surf),
      Conc.4hour = mean(Conc.4hour),
      Conc.surf.i = mean(Conc.surf.i)
    )
    
    # Hand concentration after touching surface
    Conc.h[i, ] <- Conc.h.prev * exp(-k.hand * Time.m1) - {
      param$Frac.hs * (TE.hs * Conc.h.prev * exp(-k.hand * Time.m1) -
                         TE.sh * Conc.4hour * exp(-k.surf * Time.m1))
    }
    
    
    # Surface concentration update
    Conc.s[i, ] <- Conc.4hour * exp(-k.surf * Time.m1)- {
      param$Frac.hs * (T.handarea / param$T.surfarea) *
        (TE.sh * Conc.4hour * exp(-k.surf * Time.m1) -
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
# Step 2: Model function to pack and run all the scenario in different sequence (sub-scenarios)
#==============================================================================================
run_model2_all <- function(iter) {
  
  # scenarios: Baseline + I1~I6
  scenarios <- c("Baseline", "I1", "I2", "I3", "I1+I2", "I1+I2+I3", "I2+I3")
  
  # sequence list 
  sequences <- list(
    E = c("Elevator"),
    FD = c("Frontdesk"),
    E_FD = c("Elevator", "Frontdesk"),
    TT= c("Table"),
    E_TT = c("Elevator", "Table")
  )
  
  results <- list()
  
  for (sc in scenarios) {
    for (seq_name in names(sequences)) {
      seq_val <- sequences[[seq_name]]
      out <- run_touch_sequence_2(seq_val, sc, iter)
      key <- paste0("Model2_", sc, "_", seq_name)
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
library(ggbreak)

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
  
  #Scenario order 
  full_df$Scenario <- factor(full_df$Scenario,
                             levels = c("Baseline", "I1", "I2", "I3", "I1+I2", "I2+I3", "I1+I2+I3"))
  
  return(full_df)
}

# Violin plot visualization 

plot_risk_violin <- function(df) {
  ggplot(df, aes(x = Scenario, y = Risk, fill = Scenario)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.8) +
    facet_grid(rows = vars(Sequence), scales = "free_y", switch = "x") +
    #facet_wrap(~ Sequence , scales = "fixed") +
    scale_y_continuous(
      trans = "log10",
      labels = scales::scientific
      #,breaks = c(1e-20, 1e-18, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6)
    ) +
    labs(y = "Risk (log10)", x = "Scenario") +
    theme_minimal()
}

plot_risk_violin_1 <- function(df) {
  ggplot(df, aes(x = Scenario, y = Risk, fill = Scenario)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.8) +
    scale_y_continuous(trans = "log10", labels = scales::scientific,
                       breaks = c(1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1)) +
    labs(y = "Risk (log10)", x = "Scenario") +
    theme_minimal()
}

#========================================
#Step 4: Run model 1 and plot the model 1
#========================================


#4.1: Run model
results <- run_model2_all(iter)

#Parameter saving 
debug_log_df <- do.call(rbind, debug_log_list)
write.csv(debug_log_df, "debug_output_m2_layer.csv", row.names = FALSE)


#4.2: Data prep for visualization
risk_df <- prepare_risk_plot_data(results)


#4.3: Visualization 

# Violin plot for total scenarios

windows()
plot_risk_violin(risk_df)

ggsave("Infection risk_m2_layer(all).tiff", dpi=600, dev= 'tiff', height=6, width=20, units='in')

#4.4: Intervention compare

risk_df_comp <-risk_df %>% filter(Scenario %in% c("Baseline", "I1", "I2", "I3", "I1+I2","I2+I3", "I1+I2+I3"))


# Set the threshold (1e-22)
threshold <- 1e-22

risk_df_comp_marked <- risk_df_comp %>%
  mutate(
    Risk_display = pmax(Risk, threshold),  # change the value below threshold to threshold value
    below_threshold = Risk < threshold     
  )

#extreme value information summary 
extreme_values<-risk_df_comp %>%
  filter(Risk<threshold) %>%
  group_by(Sequence, Scenario) %>%
  summarize (
    count = n(),
    min_risk= min(Risk),
    mean_risk=mean(Risk),
    .groups='drop'
  )

#Plot drawing 
plot_itv <- ggplot(risk_df_comp_marked, aes(x = Scenario, y = Risk_display, fill = Scenario)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.8) +
  
  # data below threshold line 
  geom_point(data = filter(risk_df_comp_marked, below_threshold), 
             aes(y = Risk_display), 
             shape = 25, size = 2, color = "red", fill = "red") +
  
  #Add extreme values 
  geom_text(data=extreme_values, 
            aes(x= Scenario, y=threshold,
                label=paste0("Mean:", scales::scientific (mean_risk, digits=1),
                             "\nMin:", scales::scientific (min_risk, digits=1))),
            vjust=1.2, hjust=0.5, size=2.5, color = "red",
            fontface = "bold")+
  
  facet_wrap(~ Sequence, scales = "fixed", nrow=1) +
  scale_y_continuous(trans = "log10", 
                     labels = scales::scientific,
                     breaks = c(1e-20, 1e-18, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6),
                     limits = c(threshold, max(risk_df_comp$Risk)),
                     expand = expansion(mult = c(0.1, 0))) +
  
  # Threshold line drawing
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(y = "Risk (log10)", x = "Scenario",
       caption = "Red triangles indicate values below 1e-22") +
  theme_minimal()

windows()
plot_itv

ggsave("Infection risk_m2_layer.tiff", dpi=600, dev= 'tiff', height=6, width=20, units='in',limitsize = FALSE)



#==============================================================
#Step 5: Summary Statistic extraction & save result CSV. files
#==============================================================


#Step 5.0: Data structure confirm 
str(results[["Model2_Baseline_E"]]$Conc.h)
df_test<- as.data.frame(results[["Model2_Baseline_E"]]$Conc.h)

#Step 5.1: Extract event-level values from each result matrix and reshape for analysis
prepare_event_level_data <- function (results, variable_name=c("Conc.h", "Conc.s", "Dose", "Risk")){
  df_all<-list()
  
  for (var in variable_name){
    temp_list<-list()
    for (name in names(results)) {
      res<-results[[name]]
      scenario_seq <-strsplit(name, "_")[[1]]
      scenario <-scenario_seq[2]
      sequence <-paste (scenario_seq[3:length(scenario_seq)], collapse="_")
      mat <-res[[var]] #e.g.) pull out the matrix result such as "Conc.h"
      
      df <- data.frame(
        value =as.vector(t(mat)), # Switch row and column of result matrix and make it vector in one column (read by row/ event - by - iteration ) 
        event = rep(rownames(mat), each= ncol(mat)), # put event name on each value
        scenario = scenario,
        sequence = sequence
      )
      temp_list[[name]] <-df
    }
    df_all[[var]] <-do.call(rbind, temp_list)
  }
  return(df_all)
}


#Step 5.2: Summary statistic extraction function with percent risk reduction for Risk matrix

generate_summary_stats <- function (df, value_name){
  stats <- df %>%
    group_by(scenario, sequence, event) %>%
    summarize (
      mean = mean(value, na.rm=TRUE),
      sd = sd(value, na.rm=TRUE),
      min = min(value, na.rm=TRUE),
      max = max(value, na.rm=TRUE),
      median = median(value, na.rm=TRUE),
      .groups= "drop"
    ) %>%
    mutate(variable=value_name) %>%
    select(scenario, sequence, event, variable, everything())
  
  return(stats)
}


# Step 5.3: Use summary stat function and write csv. file

all_data <- prepare_event_level_data (results)

summary_Conc_h<-generate_summary_stats(all_data$Conc.h, "Conc.h")
summary_Conc_s<-generate_summary_stats(all_data$Conc.s, "Conc.s")
summary_Dose<-generate_summary_stats(all_data$Dose, "Dose")
summary_Risk<-generate_summary_stats(all_data$Risk, "Risk")

write.csv(summary_Conc_h, "summary_Conc_h_m2.csv", row.names = FALSE)
write.csv(summary_Conc_s, "summary_Conc_s_m2.csv", row.names = FALSE)
write.csv(summary_Dose,   "summary_Dose_m2.csv", row.names = FALSE)
write.csv(summary_Risk,   "summary_Risk_m2.csv", row.names = FALSE)
