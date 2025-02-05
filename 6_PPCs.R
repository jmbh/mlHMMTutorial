# jonashaslbeck@protonmail.com; Feb 5th, 2025

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Here we perform posterior predictive checks on the 
# distributions of within-person means and SDs for all variables
# However, we only show them for the variables Happy and Sad, 
# since they are representative for all eight positive/negative
# emotions

# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------

library(mHMMbayes)
library(xtable)
library(plyr)
library(dplyr)
library(RColorBrewer)

source("0_Helpers.R")

# --------------------------------------------------------
# ---------- Load Data -----------------------------------
# --------------------------------------------------------

emotion_mHMM <- readRDS("Data/data_Rowland2020_wN_DF.RDS")


# --------------------------------------------------------
# ---------- Load Model Object ---------------------------
# --------------------------------------------------------

l_Models <- readRDS("Files/RowlandHMMs_rep1.RDS")


# --------------------------------------------------------
# ---------- Compute Empirical Within-person Means/SDs ---
# --------------------------------------------------------

# Means
pw_means_emp <- ddply(emotion_mHMM, .(subj_id), function(x) colMeans(x[,2:9], na.rm = TRUE))
# SDs
pw_sds_emp <- ddply(emotion_mHMM, .(subj_id), function(x) apply(x[,2:9], 2, function(y) sd(y, na.rm=TRUE)))

# --------------------------------------------------------
# ---------- Simulate Data & Get Predicted Means / SDs ---
# --------------------------------------------------------

# We simulate data given the six models

l_Ms <- l_SDs <- list()

for(m in 1:6) {
  
  # Obtaining parameters
  model_m <- l_Models[[m]]
  ## ----- A) Getting Fixed effects -----
  # Transition Probabilities
  if(m==1) {
    gamma_group <- matrix(1, 1, 1)
  } else {
    gamma_group <- obtain_gamma(model_m)
  }
  # Emission Distributions
  emiss_group <- obtain_emiss(model_m)
  
  ## ------ B.2) Getting Random Effects VAR from explicit parameters ------
  ## Transition
  if(m>1) {
    m_trans_vars <- apply(model_m$gamma_V_int_bar, 2, mean, na.rm=TRUE)
    m_trans_vars # Larger than the ones from the REs (see above), as Emmeke predicted
  } else {
    m_trans_vars <- matrix(1, 1, 1)
  }
  ## Emission
  # This takes the average over posterior samples of explicitly modeled variances
  m_emiss_vars <- t(sapply((model_m$emiss_varmu_bar), apply, 2, mean, na.rm = TRUE)) # as in Practical of Emmeke
  m_emiss_vars <- matrix(m_emiss_vars, nrow=8, ncol=m, byrow = FALSE)
  # Get this in the right format
  l_emiss_vars <- list()
  for(i in 1:8) l_emiss_vars[[i]] <- matrix(m_emiss_vars[i, ], )
  # Renormalize gamma, to avoid floating point issues
  gamma_group_renorm <- gamma_group / rowSums(gamma_group)
  gamma_group_renorm <- matrix(gamma_group_renorm, m, m)
  
  # Generate data
  data_m <- sim_mHMM(n_t = 240, #
                     # n = n_subj,
                     n = 1000, # We do more here than in the data, to approximate the expectation
                     data_distr = "continuous",
                     gen = list(m=m, n_dep=n_dep),
                     gamma = gamma_group_renorm, # fixed effects: matrix; we round here to avoid some floating point issue in input checks on transition matrix
                     emiss_distr = emiss_group, # fixed effects: emissions
                     var_gamma = m_trans_vars, # requires: m x m-1 matrix
                     var_emiss = l_emiss_vars) # Correct
  
  # Convert to dataframe
  data_f <- data.frame(data_m$obs)
  labels <- colnames(emotion_mHMM)
  colnames(data_f) <- labels
  
  # Compute person-wise means
  l_Ms[[m]] <- ddply(data_f, .(subj_id), function(x) colMeans(x[,2:9], na.rm = TRUE))
  l_SDs[[m]] <- ddply(data_f, .(subj_id),function(x) apply(x[,2:9], 2, function(y) sd(y, na.rm=TRUE)))
  
  # Progress
  print(m)
  
} # end for: m

# ----- Save -----
# saveRDS(l_Ms, "Files/PPC_SimMeans.RDS")
# saveRDS(l_SDs, "Files/PPC_SimSDs.RDS")

# Load
l_Ms <- readRDS("Files/PPC_SimMeans.RDS")
l_SDs <- readRDS("Files/PPC_SimSDs.RDS")


# --------------------------------------------------------
# ---------- Make Figure: Cmb Mean & SD for Happy/Sad ----
# --------------------------------------------------------

# This is Figure 9 in the paper

labelsU <- c("Happy", "Excited", "Relaxed", "Satisfied",
             "Angry", "Anxious", "Depressed", "Sad")

m_max <- 4
cols <- brewer.pal(m_max, "Set1")


sc <- 0.9
pdf("Figures/PPC_MeanSD_HappySad.pdf", width=7*sc, height = 5*sc)

# Layout
lmat <- matrix(5:8, 2, 2, byrow = TRUE)
lmat <- rbind(1:2, lmat)
lmat <- cbind(c(0, 3:4), lmat)
wh <- c(0.15, 1, 1)
lo <- layout(lmat, widths = wh, heights = wh)
# layout.show(lo)

# Labels
plotLabel("Happy", cex=1.6)
plotLabel("Sad", cex=1.6)
plotLabel("Within-person Means", srt=90, cex=1.6)
plotLabel("Within-person SDs", srt=90, cex=1.6)

# --- Means ---
for(i in c(1,8)) {
  # Get densities
  den_emp <- density(pw_means_emp[, 1+i])
  
  # Plot canvas
  par(mar=c(3,4,1,1))
  plot.new()
  plot.window(xlim=c(0, 100), ylim=c(0, 0.052))
  axis(1)
  axis(2, las=2)
  title(ylab="Density")

  if(i==1) legend("topleft",
                  legend=c("Empirical", paste0("M = ", 1:m_max)),
                  text.col = c("black", cols), bty="n", cex=1.2)
  
  # Plot Data
  lines(den_emp, lwd=2)
  # Fetch simulated w-p means
  for(m in 1:m_max) {
    pw_means_sim <- l_Ms[[m]][, -1][, i]
    den_sim <- density(pw_means_sim)
    lines(den_sim, col=cols[m], lwd=2)
  } # end: m
}


# --- SDs ---
for(i in c(1,8)) {
  # Get densities
  den_emp <- density(pw_sds_emp[, i+1])
  
  # Plot canvas
  plot.new()
  plot.window(xlim=c(0, 60), ylim=c(0, 0.4))
  axis(1)
  axis(2, las=2)
  title(ylab="Density")
  # Plot Data
  lines(den_emp, lwd=2)
  # Fetch simulated w-p means
  for(m in 1:m_max) {
    pw_SDs_sim <- l_SDs[[m]][, -1][, i]
    den_sim <- density(pw_SDs_sim)
    lines(den_sim, col=cols[m], lwd=2)
  } # end: m
  
}

dev.off()


# --------------------------------------------------------
# ---------- Plot PPC Means (all variables) --------------
# --------------------------------------------------------

# This is shown in the appendix

m_max <- 4
cols <- brewer.pal(m_max, "Set1")

# ----- For m=1:4 (for paper) -----
sc <- 1.2
pdf(paste0("Figures/PPC_pw_Means_m=1to", m_max , ".pdf"), width = 8*sc, height=4*sc)

par(mfrow=c(2, 4))

for(i in 1:8) {
  
  # Get densities
  den_emp <- density(pw_means_emp[, i+1])
  
  # Plot canvas
  par(mar=c(4,4,2,1))
  plot.new()
  plot.window(xlim=c(0, 100), ylim=c(0, 0.052))
  axis(1)
  axis(2, las=2)
  title(ylab="Density")
  title(xlab="Within-person Mean", line=2.5)
  title(main=labelsU[i], font.main=1)
  if(i==2) {
    legend("topright",
           legend=c("Empirical", paste0("M = ", 1:m_max)),
           text.col = c("black", cols), bty="n", cex=1.2)
  }
  
  # Plot Data
  lines(den_emp, lwd=2)
  # Fetch simulated w-p means
  for(m in 1:m_max) {
    pw_means_sim <- l_Ms[[m]][, -1][, i]
    den_sim <- density(pw_means_sim)
    lines(den_sim, col=cols[m], lwd=2)
  } # end: m
} # end: i

dev.off()


# --------------------------------------------------------
# ---------- Plot PPC SDs (all variables) ----------------
# --------------------------------------------------------

# This is shown in the appendix


m_max <- 4
cols <- brewer.pal(m_max, "Set1")

# ----- For m=1:4 (for paper) -----
sc <- 1.2
pdf(paste0("Figures/PPC_pw_SDs_m=1to", m_max , ".pdf"), width = 8*sc, height=4*sc)

par(mfrow=c(2, 4))

for(i in 1:8) {
  
  # Get densities
  den_emp <- density(pw_sds_emp[, i+1])
  
  # Plot canvas
  par(mar=c(4,4,2,1))
  
  plot.new()
  plot.window(xlim=c(0, 60), ylim=c(0, 0.2))
  axis(1)
  axis(2, las=2)
  title(ylab="Density")
  title(xlab="Within-person SD", line=2.5)
  title(main=labelsU[i], font.main=1)
  if(i==2) {
    legend("topright",
           legend=c("Empirical", paste0("M = ", 1:m_max)),
           text.col = c("black", cols), bty="n", cex=1.2)
  }
  
  # Plot Data
  lines(den_emp, lwd=2)
  # Fetch simulated w-p means
  for(m in 1:m_max) {
    pw_SDs_sim <- l_SDs[[m]][, -1][, i]
    den_sim <- density(pw_SDs_sim)
    lines(den_sim, col=cols[m], lwd=2)
  } # end: m
} # end: i

dev.off()







