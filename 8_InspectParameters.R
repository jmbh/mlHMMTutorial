# jonashaslbeck@protonmail.com; Feb 5th, 2025

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Here we extract parameter estimates of fixed and random
# effects estimates of the emission distribution and the
# transition probability matrix and plot associated figures
# shown in the paper and appendix

# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------

library(mHMMbayes)
library(xtable)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(scales)

source("0_Helpers.R")

# --------------------------------------------------------
# ---------- Load Model Object ---------------------------
# --------------------------------------------------------

l_Models <- readRDS("Files/RowlandHMMs_rep1.RDS")


# --------------------------------------------------------
# ---------- Extract & Plot Emission Estimates -----------
# --------------------------------------------------------

# Get plots separatelt for Means and SDs and for models with
# states 2,3,4

# This uses the function PlotEmiss from 0_Helpers.R
v_pars <- c("Means", "SDs")
for(m in 2:4) {
  for(par in 1:2) {
    PlotEmiss(object = l_Models,
              m = m,
              parameter=v_pars[par])
  }
}

# --------------------------------------------------------
# ---------- Extract & Plot Transition Estimates ---------
# --------------------------------------------------------

# Show them for states 2,3,4

for(m in 2:4) {
  
  n_dep <- 8
  N <- 125
  model_m <- l_Models[[m]]
  
  # ------ Getting the Parameters ------
  ## Fixed Effects
  gamma_group <- obtain_gamma(model_m)
  ## Random Effects
  gamma_subj <- obtain_gamma(model_m, level = "subject")
  a_gamma <- array(NA, dim=c(m, m, N))
  for(i in 1:N) a_gamma[, , i] <- gamma_subj[[i]]
  
  # ------ Plotting ------
  pdf(paste0("Figures/Transition_m=", m, ".pdf"), width=8.5, height=4)
  
  # Layout
  lo <- layout(matrix(1:2, ncol=2), widths = c(.4, .6))
  # layout.show(lo)
  
  ### Left Panel: Transition Matrix
  plt_h <- plotHeat(gamma_group, cex.axis=1.2, cex.ft=1.2)
  
  ### Right Panel: Barplot with Random effects
  # Generate x_labels
  v_labels <- rep(NA, m^2)
  cnt <- 1
  for(m1 in 1:m) for(m2 in 1:m) {
    v_labels[cnt] <- paste0("S", m1, " to ", "S", m2)
    cnt <- cnt+1
  }
  
  # Barplot
  par(mar=c(4.5,3,3,1))
  bp <- barplot(as.numeric(t(gamma_group)),
                ylim=c(0, 1),
                col = as.character(plt_h[m:1,]),
                # names.arg = v_labels,
                cex.names=0.7, axes=FALSE)
  axis(1, v_labels, las=2, at=bp)
  axis(2, las=2)
  
  # Add Random Effects
  N_disp <- 50
  cnt <- 1
  for(m1 in 1:m) {
    for(m2 in 1:m) {
      points(rep(bp[cnt], N_disp), a_gamma[m1, m2, 1:N_disp],
             pch=21, col="black",
             bg=as.character(plt_h[m:1,])[cnt], cex=0.75)
      cnt <- cnt + 1
    }
  }
  
  dev.off()
  
} # end loop: m states








