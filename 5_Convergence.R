# jonashaslbeck@protonmail.com; Feb 5th, 2025

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Here we inspect the log-likelihood (model fit) and the
# AIC for all considered models


# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------

library(mHMMbayes)
library(xtable)
library(RColorBrewer)
library(coda) 

source("0_Helpers.R")


# --------------------------------------------------------
# ---------- Load Model Object for each Chain ------------
# --------------------------------------------------------

Nchain <- 4 # Number of chains

l_Models <- list()
for(chain in 1:Nchain) {
  l_Models[[chain]] <- readRDS(paste0("Files/RowlandHMMs_rep", chain, ".RDS"))
}

# --------------------------------------------------------
# ---------- Trace Plots ---------------------------------
# --------------------------------------------------------
# Just show two: one for fixed transition matrix, one for fixed

# These are shown in the appendix of the paper

# Select Colors
cols <- brewer.pal(4, "Set1")

sc <- 1.3
pdf(paste0("Figures/Traceplots_CMB_4chains.pdf"), width=4*sc, height=6.5*sc)

# Make Layout
lmat <- matrix(9:20, 6, 2, byrow = TRUE)
lmat <- rbind(7:8, lmat)
lmat <- cbind(c(0, 1:6), lmat)
lo <- layout(mat=lmat, widths = c(.25, 1, 1), heights = c(.2, rep(1, 6)))
# Fill in Labels
for(m in 1:6) plotLabel(paste0("M = ", m))
plotLabel(expression("Logit: " ~ S[1] ~ " to " ~ S[2]))
plotLabel(expression("Emission Mean 'Happy' in "~ S[max]))

# Loop in Data
par(mar=c(2.2,2,1,1))
for(m in 1:6) {
  
  # Trace for Transition Matrix
  if(m>1) {
    plot.new()
    plot.window(xlim=c(1,N), ylim=c(-4,0.5))
    axis(1)
    axis(2, las=2)
    for(chain in 1:Nchain) {
      lines(l_Models[[chain]][[m]]$gamma_int_bar[,1], col=cols[chain])
    }
  } else {
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    legend("center", legend = paste0("Chain", 1:4),
           text.col = cols, bty = "n", cex = 1.2)
    
  }
  # Trace for Emission means
  plot.new()
  plot.window(xlim=c(1,N), ylim=c(40,80))
  axis(1)
  axis(2, las=2)
  for(chain in 1:Nchain) {
    lines(l_Models[[chain]][[m]]$emiss_mu_bar$happy[, m], col=cols[chain])
  }
  
} # end for: m

dev.off()


# --------------------------------------------------------
# ---------- Compute Gelman/Ruben (Rhat) Statistics ------
# --------------------------------------------------------

# ----- Compute GR Statistic -----
n_dep <- 8
l_RG <- list()

for(m in 1:6) {
  model_list <- list()
  for(chain in 1:Nchain) model_list[[chain]] <- l_Models[[chain]][[m]]
  l_RG[[m]] <- f_GR(model_list) # from: 0_Helpers.R
} # end for: m

# ----- Save -----
saveRDS(l_RG, file="Files/Convergence_GMStat_4chains.RDS")

# ----- Summarize -----
# Get means per m and parameter type
RG_agg <- lapply(l_RG, function(x) round(c(mean(x[[1]]), mean(x[[2]])),2) )
m_RG_agg <- do.call(cbind, RG_agg)
colnames(m_RG_agg) <- paste0("M = ", 1:6)
rownames(m_RG_agg) <- c("Transition", "Emission")

xtab <- xtable(m_RG_agg, digits = rep(2, 7))
# Save LateX Table into folder "Files"
print(xtab, type = "latex", file = "Files/table_RG.tex", include.rownames = TRUE)





