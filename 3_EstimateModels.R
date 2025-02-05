# jonashaslbeck@protonmail.com; Feb 5th, 2025

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# We specify priors and starting values and fit mlHMM with
# a sequence of m=1:6 states


# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------

library(mHMMbayes)


# --------------------------------------------------------
# ---------- Load Data -----------------------------------
# --------------------------------------------------------

emotion_mHMM <- readRDS("Data/data_Rowland2020_wN_DF.RDS")


# --------------------------------------------------------
# ---------- Estimate Sequence ---------------------------
# --------------------------------------------------------

n_dep <- 8 # Number of dependent variables
J <- 2000 # Number of Gibbs samples
burn_in <- 500 # Burnin for Gibbs

# Specify m-sequence
m <- 1:6

# ---------- Loop over m ----------

# We run four separate chains here, in order to be able to assess convergence

for(chain in 1:4) {
  
  # Storage
  l_Models <- list()
  
  # Reproducibility
  set.seed(chain)
  
  for(m in 1:6) {
    
    # ----- Specify Starting Values ------
    
    ### Starting values: Transition probabilities
    start_gamma <- matrix(c(1), byrow = TRUE, ncol = m, nrow = m)
    
    ### Starting values: Emission distribution
    start_emiss <- list()
    # Positive Valence
    i_Mean_PV <- seq(0, 100, length=m+2)[-c(1, m+2)] # Means
    i_SD_PV <- rep(15, m) # SDs
    for(i in 1:4) start_emiss[[i]] <- cbind(i_Mean_PV, i_SD_PV)
    # Negative Valence
    i_Mean_NV <- seq(0, 50, length=m+2)[-c(1, m+2)]     # Means
    i_SD_NV <- rep(10, m)
    for(i in 5:8) start_emiss[[i]] <- cbind(i_Mean_NV, i_SD_NV)
    
    # ----- Specify Priors ------
    
    ## Priors: specifying weakly informative prior for continuous emission distributions
    emotion_prior_emiss <- prior_emiss_cont(
      gen = list(m = m, n_dep = n_dep),
      emiss_mu0 = list(matrix(i_Mean_PV, nrow = 1),  # happy; this is the prior in BOTH states
                       matrix(i_Mean_PV, nrow = 1),  # excited
                       matrix(i_Mean_PV, nrow = 1),  # relaxed
                       matrix(i_Mean_PV, nrow = 1),  # satisfied
                       matrix(i_Mean_NV, nrow = 1),  # angry
                       matrix(i_Mean_NV, nrow = 1),  # anxious
                       matrix(i_Mean_NV, nrow = 1),  # depressed
                       matrix(i_Mean_NV, nrow = 1)), # sad
      emiss_K0 = rep(list(1), n_dep), # number of people the prior is based on
      emiss_V = rep(list(rep(5^2, m)), n_dep), # prior between subj variance
      # emiss_V = rep(list(rep(10^2, m)), n_dep), # prior between subj variance
      emiss_nu = rep(list(1), n_dep), # DF between subj variance
      emiss_a0 = rep(list(rep(1.5, m)), n_dep), # shape: emission
      emiss_b0 = rep(list(rep(20, m)), n_dep), # scale: emission
    )
    
    # ----- Estimate ml-HMM ------
    
    l_Models[[m]] <- mHMM(s_data = emotion_mHMM,
                          data_distr = "continuous",
                          gen = list(m = m, n_dep = n_dep),
                          start_val = c(list(start_gamma), start_emiss),
                          emiss_hyp_prior = emotion_prior_emiss,
                          mcmc = list(J = J, burn_in = burn_in))
    
    # Print Progress
    prog <- paste0("Rep=", chain, "; Model m=", m)
    print(prog)
    
  } # end for m
  
  # ----- Save ------
  saveRDS(l_Models, paste0("Files/RowlandHMMs_rep", chain, ".RDS"))
  
} # end for Repetition



