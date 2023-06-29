# A Computational Model of Depression

# This code produces a computational model of depression. Section 1 sets up and
# builds the model step by step. First, the differential equations and initial
# and default values are defined. Then, the function that carries out the 
# simulation is defined. Finally, an example of how the model can be run and
# visualized is provided.

# Section 2 contains the code for running and visualizing the simulations that
# were reported. Simulation 1a simulates the system for one year (N=10,000)
# and checks the one-year prevalence rate. Simulation repeats simulation 1a
# 10 times, checking the prevalence rate for each simulated dataset. Simulation
# 2 simulates the system for two years (N=10,000) and checks the average length
# of a depressive episode. Simulation 3 picks a subset from Simulation 2 and 
# plots mood against environmental stimuli.


# SECTION 1=====================================================================


# Differential Equations--------------------------------------------------------

# Mood
dM_dt <- function(M,
                  B,
                  E,
                  k_B_M,
                  h_B_M,
                  r_M) {
  r_M * (1/(1 + exp(- k_B_M * (B + h_B_M + E))) - M) + rnorm(1, mean = 0, sd = 10)
}


# Behavioral Activity
dB_dt <- function(B,
                  M,
                  k_M_B,
                  h_M_B,
                  r_B) {
  r_B * (1/(1 + exp(- k_M_B * (M + h_M_B))) - B) + rnorm(1, mean = 0, sd = 10)
}


# Cognitive Style
dC_dt <- function(C,
                  weightedM,
                  cr_M_C,
                  r_M_C) {
  ifelse(weightedM >= cr_M_C,
         abs(r_M_C * weightedM),
         ifelse(weightedM <= -cr_M_C,
                -abs(r_M_C * weightedM), 0))
}


# Defaults----------------------------------------------------------------------

# Parameter defaults
pars_default <- list("B" = list("r_B" = 0.05, # intrinsic rate of B
                                "k_M_B" = 20, # threshold at which M leads to B
                                "h_M_B" = 0.4), # rate at which M leads to B
                     "M" = list("r_M" = 0.05, # intrinsic rate of M
                                "k_B_M" = 20, # threshold at which B leads to M
                                "h_B_M" = 0.4), # rate at which B leads to M
                     "C" = list("cr_M_C" = 0.3,
                                "r_M_C" = 0.5))


# Initial defaults
initial_default <- list("B" = 0,
                        "M" = 0,
                        "C" = 0,
                        "E" = 0)

# Simulation defaults
sim_default <- list("stepsize" = .001)




# Simulate----------------------------------------------------------------------
simDepression <- function(time, # integer vector 1:n, indicating the time interval, where 1 is one "day"
                          stepsize = NULL, # stepsize; <= 1
                          parameters = NULL, # specify parameter values for the model
                          initial = NULL, # specify initial values for state variables
                          pbar = TRUE) # progress bar
  {
  
  # Setup Step 0: Overwrite default with specified parameters/initial values
  
  ## Overwrite default parameters, if specified
  PS <- pars_default
  
  # Parameters specified?
  if (!is.null(parameters)) {
    
    names_list_def <- names(pars_default)
    names_list_spec <- names(parameters)
    n_spec <- length(names_list_spec)
    
    # Loop over specified upper level list entries
    for (i in 1:n_spec) {
      pars_spec_i <- parameters[[names_list_spec[i]]]
      n_pars_spec_i <- length(pars_spec_i)
      names_list_spec_j <- names(pars_spec_i)
      for(j in 1:n_pars_spec_i) {
        PS[[names_list_spec[i]]][[names_list_spec_j[j]]] <- parameters[[names_list_spec[i]]][[names_list_spec_j[j]]]
      }
    }
  }
  
  ## Overwrite default starting values, if specified
  INI <- initial_default
  
  # Initial values specified?
  if (!is.null(initial)) {
    
    names_list_ini_def <- names(initial_default)
    names_list_ini_def <- names(initial)
    n_spec_ini <- length(names_list_ini_def)
    
    for (i in 1:n_spec_ini) {
      INI[[names_list_ini_def[i]]] <- initial[[names_list_ini_def[i]]]
    }
    
  }
  
  ## Specify simulation parameters, if not previously specified
  if (is.null(stepsize)) {
    stepsize <- sim_default$stepsize
  }
  
  # Setup Step 1: Import time scale and default parameters
  
  range_time <- range(time) # range of time steps, 1:number of specified iterations
  
  
  # Setup Step 2: Specify components
  
  # Set initial values
  B <- INI$B # behavioral activity
  M <- INI$M # mood
  C <- INI$C # cognitive style
  E <- INI$E # environment
  
  
  # Setup Step 3: Specify additional parameters
  
  # Initial Cognitive Style parameters 
  PS$B$k_M_B <- 20 - 10 * 0.10^C
  PS$B$h_M_B <- 0.25^C
  
  # Create storage
  outmat <- matrix(NA, nrow = length(time), ncol = 4,
                   dimnames = list(NULL, c("B", "M", "C", "E")))
  
  # Save initial values
  outmat[1, ] <- c(B, M, C, E)
  weightedM <- 0
  
  # Track time with a time-tracker
  obs_tracker <- 2 # the first observation has already been taken, so next save 
  # point is the second one
  
  # Simulation Step 1: 'For Loop' that carries out simulation
  
  # Setup progress bar
  if (pbar == TRUE) {
    pb <- txtProgressBar(min = 1, max = length(time), initial = 1,
                         char = "-", style = 3)
  }
  
  # The For Loop
  for (obs_tracker in 2:length(time)) {
    
    # Simulation Step 1: Fast Model Situations
    
    # Update other variables based on past
    
    ## Create new values
    
    # Behavioral activity
    Bnew <- B + dB_dt(B = B,
                      M = M,
                      k_M_B = PS$B$k_M_B,
                      h_M_B = PS$B$h_M_B,
                      r_B = PS$B$r_B) * stepsize
    
    # Mood 
    Mnew <- M + dM_dt(M = M,
                      B = B,
                      E = E,
                      k_B_M = PS$M$k_B_M,
                      h_B_M = PS$M$h_B_M,
                      r_M = PS$M$r_M) * stepsize
    
    # Environment
    # The probability of E switching increases the longer it has stayed the same
    if (obs_tracker > 1) {
      prev_E <- outmat[obs_tracker - 1, "E"]
      current_chunk_length <- 0
      for (i in (obs_tracker - 1):1) {
        if (outmat[i, "E"] == prev_E) {
          current_chunk_length <- current_chunk_length + 1
        } else {
          break
        }
      }
      if (current_chunk_length >= 7) {
        switch_prob <- 0.1  # High probability of switching if current chunk is long
      } else {
        switch_prob <- 0.1 * (7 - current_chunk_length)  # Decrease probability as chunk gets longer
      }
      if (runif(1) < switch_prob) {
        Enew <- sample(c(-1, 0, 1), size = 1, prob = c(1/3, 1/3, 1/3), replace = TRUE)
      } else {
        Enew <- prev_E
      }
    } else {
      Enew <- E
    }
    
    
    # Overwrite current values
    B <- Bnew
    M <- Mnew
    E <- Enew
    
    # Update weightedM
    
    weights <- sort((1/(1:30))/sum(1/(1:30)), decreasing = FALSE) # create
    # weights for determining the relative importance of the past month's mood,
    # the weights follow y = 1/x
    
    if (i >= 30) {
      weightedM <- weighted.mean(x = outmat[(i-29):i, "M"],
                                 weights = weights)
    } else {
      weightedM <- 0
    }
    
    # Update cognitive style
    Cnew <- dC_dt(C = C,
                  weightedM = weightedM,
                  cr_M_C = PS$C$cr_M_C,
                  r_M_C = PS$C$r_M_C) # + C
    
    # Overrride cognitive style
    C <- Cnew
    
    # Update parameters determined by cognitive style
    PS$B$k_B_M <- 20 - 10 * 0.1^C
    PS$B$h_B_M <- 1.1^C
    
    
    # Save current values to the output matrix
    outmat[obs_tracker, ] <- c(B, M, C, E)
    
    # Update observation tracker
    obs_tracker <- obs_tracker + 1
    
    # Update progress bar
    if (pbar == TRUE) {
      setTxtProgressBar(pb, obs_tracker)
    }
    
  }
  
  # Simulation Step 2: Specify function's output
  
  outmat <- as.data.frame(outmat)
  
  outlist <- list("outmat" = outmat,
                  "input" = list("parameters" = PS,
                                 "initial" = INI))
  return(outlist)
  
}


# Run and visualize model: Example----------------------------------------------

out <- simDepression(time = 1:2000)

# Store results
results <- out$outmat
plot.new()
plot.window(xlim = c(0, length(results$M)), ylim = c(-1, 1))
axis(1); axis(2, las = 2)
mtext("Days", side = 1, line = 3, cex = 1.25)
lines(results$M, col = "#E56399", lwd = 1.5)
lines(results$B, col = "#7F96FF", lwd = 1.5)
lines(results$C, col = "#320E3B", lwd = 1.5)
legend("topright", legend = c("Mood", "Behavior", "Cognitive Style"), 
       bty = "n", cex = 1, fill = c("#E56399", "#7F96FF", "#320E3B"))

# ==============================================================================







# SECTION 2=====================================================================


# Simulation 1a-----------------------------------------------------------------

# Simulate the system for a year for N=10,000 (vary starting values for M, B, C)
# track who was above/below a certain threshold for two weeks or longer

# Modify the initial_default list with mean values
initial_default <- list("B" = list("mean" = 0, "sd" = 0.25),
                        "M" = list("mean" = 0, "sd" = 0.25),
                        "C" = list("mean" = 0, "sd" = 0.25), 
                        "E" = 0)

# Simulation function with modifications
simDepressionModified <- function(n, time, stepsize = NULL, 
                                  parameters = NULL, initial = NULL, 
                                  progbar = TRUE) {
  # Create a list to store the output matrices for each individual
  outputList <- vector("list", n)
  
  # Create a text progress bar
  if (progbar) {
    progb <- txtProgressBar(min = 0, max = n, style = 3)
  }
  
  # Loop over the number of individuals
  for (i in 1:n) {
    # Draw random initial values from normal distributions
    initialModified <- list("B" = rnorm(1, mean = initial$B$mean, 
                                        sd = initial$B$sd)[1],
                            "M" = rnorm(1, mean = initial$M$mean, 
                                        sd = initial$M$sd)[1],
                            "C" = rnorm(1, mean = initial$C$mean, 
                                        sd = initial$C$sd)[1],
                            "E" = initial$E)
    
    # Simulate depression dynamics for the current individual
    simulation <- simDepression(time = time, 
                                stepsize = stepsize, 
                                parameters = parameters, 
                                initial = initialModified, 
                                pbar = FALSE)
    
    # Store the output matrix for the current individual
    outputList[[i]] <- simulation$outmat
    
    # Update progress bar
    if (progbar == TRUE) {
      setTxtProgressBar(progb, i)
    }
    
  }
  
  return(outputList)
}

# Simulate 10000 individuals over one year
n1 <- 10000
time1 <- 1:365

# Define the initial values with mean and standard deviation
initial_modified <- list("B" = list("mean" = initial_default$B$mean, 
                                    "sd" = initial_default$B$sd),
                         "M" = list("mean" = initial_default$M$mean, 
                                    "sd" = initial_default$M$sd),
                         "C" = list("mean" = initial_default$C$mean, 
                                    "sd" = initial_default$C$sd),
                         "E" = initial_default$E)


# Simulate
output1 <- simDepressionModified(n = n1, time = time1, 
                                 initial = initial_modified)

# Set depression threshold
threshold_mood <- -0.5

# Calculate the proportion of individuals classified as depressed
depressedCount <- 0
for (i in 1:n1) {
  # Check if there is a consecutive period of 14 days or longer with mood below -0.75
  mood <- output1[[i]][, "M"] < threshold_mood
  depressedPeriods <- rle(mood)
  
  # Check if there is a period of 14 consecutive days or longer
  if (any(depressedPeriods$lengths[depressedPeriods$values] >= 14)) {
    depressedCount <- depressedCount + 1
  }
}
depressedCount / n1  # Proportion depressed in one year


# Prep for plotting
output_sub <- output1[sample(1:length(output1), size = 100, replace = FALSE)]

# Extract values from column "M" for each participant
column_M <- sapply(output_sub, function(mat) mat[, "M"])

# Create a new matrix with values from column "M"
mood_mat <- matrix(unlist(column_M), nrow = 365, 
                   ncol = 100)

# Create the depression matrix
depression_mat <- ifelse(mood_mat < -0.5, 1, 0)

# Only count episodes of 2 weeks or longer
depression_mat_dep <- depression_mat
# Iterate over each column (participant)
for (col in 1:ncol(depression_mat)) {
  # Initialize variables to track the start and end of the current episode
  episode_start <- 1
  episode_end <- 1
  
  # Iterate over each row (day)
  for (row in 2:nrow(depression_mat)) {
    # Check if the current day is part of the current episode
    if (depression_mat[row, col] == 1 && depression_mat[row - 1, col] == 1) {
      # Update the episode end
      episode_end <- row
    } else {
      # Check the length of the current episode
      episode_length <- episode_end - episode_start + 1
      
      # If the episode is shorter than 14 days, assign 0s
      if (episode_length < 14) {
        depression_mat_dep[episode_start:episode_end, col] <- 0
      }
      
      # Start a new episode
      episode_start <- row
      episode_end <- row
    }
  }
  
  # Check the length of the last episode in the column
  episode_length <- episode_end - episode_start + 1
  
  # If the last episode is shorter than 14 days, assign 0s
  if (episode_length < 14) {
    depression_mat_dep[episode_start:episode_end, col] <- 0
  }
}

# Transform
depression_mat_t <- t(depression_mat_dep)

# Plot
my_palette <- c("#FFFFFF", "#C38D94")
image(1:ncol(depression_mat_t), 1:nrow(depression_mat_t),
      t(depression_mat_t), col = my_palette,
      axes = FALSE, xlab = "", ylab = "")
tick_positions <- seq(1, ncol(depression_mat_t), length.out = 13)
image(1:ncol(depression_mat_t), 1:nrow(depression_mat_t),
      t(depression_mat_t), col = my_palette,
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = tick_positions, labels = FALSE, tck = -0.02, lwd = 0.5)
axis(3, at = tick_positions, labels = FALSE, tck = 0, lwd = 0, line = -0.2)
text(seq(15, 345, length.out = 12), -12, cex = 0.8,
     1:12, xpd = TRUE, srt = 0, adj = 0.5)
mtext("One-Year Period", side = 1, line = 2.5)
title(main = "Prevalence of depressive episodes over a one-year period (N=100)",
      cex.main = 0.9)
legend(legend = c("In a depressive episode"), fill = "#C38D94",
       bty = "n", title = "", cex = 0.8, xpd = TRUE,
       x = 125, y = 125
)




# Simulation 1b-----------------------------------------------------------------

# Simulate the system for N=10,000 10 times to extract 10 prevalence rates

n_simulations <- 10
prevalence_rates <- numeric(n_simulations)
for (sim in 1:n_simulations) {
  
  n <- 10000
  time <- 1:365
  output <- simDepressionModified(n = n, time = time, initial = initial_modified)
  
  # Calculate the proportion of individuals classified as depressed
  depressedCount <- 0
  
  for (i in 1:n) {
    # Check if there is a consecutive period of 14 days or longer with mood below -0.75
    mood <- output[[i]][, "M"] < threshold_mood
    depressedPeriods <- rle(mood)
    
    # Check if there is a period of 14 consecutive days or longer
    if (any(depressedPeriods$lengths[depressedPeriods$values] >= 14)) {
      depressedCount <- depressedCount + 1
    }
  }
  # Calculate prevalence rate for the current simulation
  prevalence_rates[sim] <- depressedCount / n
}
prevalence_rates







# Simulation 2------------------------------------------------------------------

# Simulate the system for two years N=7000
# Check the average length of depressive episodes (episodes longer than 14 days)

# Sample size + time
n2 <- 10000 
time2 <- 1:730

# Simulate
output2 <- simDepressionModified(n = n2, time = time2, 
                                 initial = initial_modified) 

# Pick the people who had a depressive episode in the simulation
# Check the average length of their depressive episode(s)

# Create a subset of dataframes for individuals with depressive episodes
output_depressed <- lapply(output2, function(df) {
  # Find consecutive days where mood is below the threshold
  consecutive_days <- rle(df$M < threshold_mood)
  
  # Get the lengths of consecutive episodes where mood is below the threshold
  episode_lengths <- consecutive_days$lengths[consecutive_days$values]
  
  # Keep only the episodes that lasted 14 days or longer
  long_episodes <- episode_lengths[episode_lengths >= 14]
  
  # Return the dataframe if there is at least one long episode
  if (length(long_episodes) > 0) {
    return(df)
  }
  
  return(NULL)
})

# Calculate the average length of depressive episodes for each individual
average_lengths <- numeric(length(output_depressed))

for (i in seq_along(output_depressed)) {
  df <- output_depressed[[i]]
  
  # Find consecutive days where mood is below the threshold
  consecutive_days <- rle(df$M < threshold_mood)
  
  # Get the lengths of consecutive episodes where mood is below the threshold
  episode_lengths <- consecutive_days$lengths[consecutive_days$values]
  
  # Keep only the episodes that lasted 14 days or longer
  long_episodes <- episode_lengths[episode_lengths >= 14]
  
  # Calculate the average length of episodes
  if (length(long_episodes) == 1) {
    average_lengths[i] <- long_episodes
  } else if (length(long_episodes) > 1) {
    average_lengths[i] <- mean(long_episodes)
  }
}
# The average length of a depressive episode for each individual
average_lengths
# The average length of a depressive episode for each depressed individual
average_lengths[which(average_lengths > 0)]
mean(average_lengths[which(average_lengths > 0)]) # "grand" mean
median(average_lengths[which(average_lengths > 0)]) # "grand" median

# Histogram
p <- ggplot(data.frame(lengths = average_lengths[which(average_lengths > 0)])) +
  geom_histogram(aes(x = lengths), bins = 30, fill = "steelblue", color = "white") +
  xlim(14, 730) +
  geom_vline(xintercept = median(filtered_lengths), linetype = "dashed", 
             color = "#FF6961", ymin = 0, ymax = max(table(filtered_lengths))) +
  labs(title = "Histogram of Average Depressive Episode Length (N=841)",
       x = "Episode Length", y = "Count") +
  # theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.background = element_blank())
p + coord_cartesian(ylim = c(0, 200),
                    expand = expansion(add = c(0, 0)))

# Prep for plotting
output_sub_2 <- output2[1:15]
depr_sub_2 <- rep(NA, length(output_sub_2))
for (i in 1:length(output_sub_2)) {
  # Check if there is a consecutive period of 14 days or longer with mood below -0.75
  mood <- output_sub_2[[i]][, "M"] < threshold_mood
  depressedPeriods <- rle(mood)
  
  # Check if there is a period of 14 consecutive days or longer
  if (any(depressedPeriods$lengths[depressedPeriods$values] >= 14)) {
    depr_sub_2[i] <- 1
    # depressedCount <- depressedCount + 1
  } else {
    depr_sub_2[i] <- 0
  }
}

# Create a df with mood values + ID + day
mood_df <- data.frame(mood = rep(NA, 730*15),
                      day = rep(1:730, 15),
                      ID = rep(1:15, each = 730))
mood_2 <- c()
for (i in 1:length(output_sub_2)) {
  mood_2 <- c(mood_2, output_sub_2[[i]][, "M"])
}
mood_df$mood <- mood_2
mood_df$depr <- rep(depr_sub_2, each = 730)

# Plot
ggplot(mood_df, aes(x = as.numeric(day), y = mood, group = ID)) +
  geom_line(aes(color = factor(depr))) +
  geom_hline(yintercept = threshold_mood, color = "darkred", linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 730, by = 100), expand = c(0, 0)) +
  labs(x = "Day", y = "Mood") +
  ggtitle("Mood fluctuations over a one-year period (N=15)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.margin = unit(0.2, "cm")
  ) +
  scale_color_manual(values = c("lightgrey", "#424242")) +
  guides(color = FALSE)








# Simulation 3------------------------------------------------------------------

# Pick individuals from previous simulation and plot their mood against
# environment values

# Plot P15
# Determine colors
p15 <- output2[15][[1]]
p15$day <- 1:nrow(p15)
cols <- rep(NA, 730)
for (i in 1:730) {
  if (p15$E[i] == -1) {
    cols[i] <- "lightgrey"
  } else if (p15$E[i] == 0) {
    cols[i] <- "white" 
  } else {
    cols[i] <- "lightyellow"
  }
}
ggplot(p15, aes(x = day, y = M)) +
  geom_rect(mapping = aes(xmin = 0:729, xmax = 1:730,
                          ymin = rep(-1, 730), ymax = rep(1, 730),
                          fill = factor(E)), color = cols) +
  geom_line(color = "black") +
  geom_hline(yintercept = threshold_mood, color = "black", 
             linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 730, by = 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("lightgrey", "white", "lightyellow"),
                    labels = c("Negative", 
                               "Neutral", 
                               "Positive")) +
  labs(x = "Day", y = "Mood", fill = "Environment") +
  ggtitle("Mood Fluctuations (N=1)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.margin = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black")))


# Plot p11
p11 <- output2[11][[1]]
p11$day <- 1:nrow(p11)
cols <- rep(NA, 730)
for (i in 1:730) {
  if (p11$E[i] == -1) {
    cols[i] <- "lightgrey"
  } else if (p11$E[i] == 0) {
    cols[i] <- "white" 
  } else {
    cols[i] <- "lightyellow"
  }
}
ggplot(p11, aes(x = day, y = M)) +
  geom_rect(mapping = aes(xmin = 0:729, xmax = 1:730,
                          ymin = rep(-1, 730), ymax = rep(1, 730),
                          fill = factor(E)), color = cols) +
  geom_line(color = "black") +
  geom_hline(yintercept = threshold_mood, color = "black", 
             linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 730, by = 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("lightgrey", "white", "lightyellow"),
                    labels = c("Negative", 
                               "Neutral", 
                               "Positive")) +
  labs(x = "Day", y = "Mood", fill = "Environment") +
  ggtitle("Mood Fluctuations (N=1)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.margin = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black")))


# Plot p1
p1 <- output2[1][[1]]
p1$day <- 1:nrow(p1)
cols <- rep(NA, 730)
for (i in 1:730) {
  if (p1$E[i] == -1) {
    cols[i] <- "lightgrey"
  } else if (p1$E[i] == 0) {
    cols[i] <- "white" 
  } else {
    cols[i] <- "lightyellow"
  }
}
ggplot(p1, aes(x = day, y = M)) +
  geom_rect(mapping = aes(xmin = 0:729, xmax = 1:730,
                          ymin = rep(-1, 730), ymax = rep(1, 730),
                          fill = factor(E)), color = cols) +
  geom_line(color = "black") +
  geom_hline(yintercept = threshold_mood, color = "black", 
             linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 730, by = 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("lightgrey", "white", "lightyellow"),
                    labels = c("Negative", 
                               "Neutral", 
                               "Positive")) +
  labs(x = "Day", y = "Mood", fill = "Environment") +
  ggtitle("Mood Fluctuations (N=1)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.margin = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black")))

# Plot p1
p1 <- output2[1][[1]]
p1$day <- 1:nrow(p1)
cols <- rep(NA, 730)
for (i in 1:730) {
  if (p1$E[i] == -1) {
    cols[i] <- "lightgrey"
  } else if (p1$E[i] == 0) {
    cols[i] <- "white" 
  } else {
    cols[i] <- "lightyellow"
  }
}
ggplot(p1, aes(x = day, y = M)) +
  geom_rect(mapping = aes(xmin = 0:729, xmax = 1:730,
                          ymin = rep(-1, 730), ymax = rep(1, 730),
                          fill = factor(E)), color = cols) +
  geom_line(color = "black") +
  geom_hline(yintercept = threshold_mood, color = "black", 
             linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 730, by = 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("lightgrey", "white", "lightyellow"),
                    labels = c("Negative", 
                               "Neutral", 
                               "Positive")) +
  labs(x = "Day", y = "Mood", fill = "Environment") +
  ggtitle("Mood Fluctuations (N=1)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.margin = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black")))



# Plot p1
p1 <- output2[1][[1]]
p1$day <- 1:nrow(p1)
cols <- rep(NA, 730)
for (i in 1:730) {
  if (p1$E[i] == -1) {
    cols[i] <- "lightgrey"
  } else if (p1$E[i] == 0) {
    cols[i] <- "white" 
  } else {
    cols[i] <- "lightyellow"
  }
}
ggplot(p1, aes(x = day, y = M)) +
  geom_rect(mapping = aes(xmin = 0:729, xmax = 1:730,
                          ymin = rep(-1, 730), ymax = rep(1, 730),
                          fill = factor(E)), color = cols) +
  geom_line(color = "black") +
  geom_hline(yintercept = threshold_mood, color = "black", 
             linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 730, by = 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("lightgrey", "white", "lightyellow"),
                    labels = c("Negative", 
                               "Neutral", 
                               "Positive")) +
  labs(x = "Day", y = "Mood", fill = "Environment") +
  ggtitle("Mood Fluctuations (N=1)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.margin = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black")))



# Plot p4
p4 <- output2[4][[1]]
p4$day <- 1:nrow(p4)
cols <- rep(NA, 730)
for (i in 1:730) {
  if (p4$E[i] == -1) {
    cols[i] <- "lightgrey"
  } else if (p4$E[i] == 0) {
    cols[i] <- "white" 
  } else {
    cols[i] <- "lightyellow"
  }
}
ggplot(p4, aes(x = day, y = M)) +
  geom_rect(mapping = aes(xmin = 0:729, xmax = 1:730,
                          ymin = rep(-1, 730), ymax = rep(1, 730),
                          fill = factor(E)), color = cols) +
  geom_line(color = "black") +
  geom_hline(yintercept = threshold_mood, color = "black", 
             linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 730, by = 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("lightgrey", "white", "lightyellow"),
                    labels = c("Negative", 
                               "Neutral", 
                               "Positive")) +
  labs(x = "Day", y = "Mood", fill = "Environment") +
  ggtitle("Mood Fluctuations (N=1)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.margin = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black")))



