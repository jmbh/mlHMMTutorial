# jonashaslbeck@protonmail.com; May 14th, 2025

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
  
  # ------ Plotting ------
  pdf(paste0("Figures/Transition_m=", m, ".pdf"), width=8.5, height=4)
  
  PlotTrans(object = l_Models,
            m = m,
            N_disp = 50) # number of random effects shown in the figure
  
  
  dev.off()
  
} # end loop: m states


# ----- Descriptives reported in paper -----
# Ranges in random effects
round(apply(a_gamma, 1:2, min), 2)
round(apply(a_gamma, 1:2, max), 2)

# Compute stationary distribution
P <- a_gamma[, , 1]
for(i in 1:3) P[i,] <- P[i,]/sum(P[i,])

# ---------------------------
# --- Compute Numerically ---
# ---------------------------

# Start with a uniform distribution
nIter <- 20
m_pi <- matrix(NA, nIter, 3)
m_pi[1, ] <- c(1/3, 1/3, 1/3)

# Iterate until convergence
for (i in 1:(nIter-1)) {
  m_pi[i+1, ] <- m_pi[i, ] %*% P
  if (max(abs(m_pi[i+1, ] - m_pi[i, ])) < 1e-10) break
  # print(i)
}
m_pi[i, ]

plot(m_pi[, 2])

# ---------------------------
# --- Compute Analytically ---
# ---------------------------

getStat <- function(P) {
  # Transpose P because eigen() finds right eigenvectors
  eig <- eigen(t(P))
  # Find the eigenvector corresponding to eigenvalue 1
  i <- which(abs(eig$values - 1) < 1e-3)  # tolerance for numerical error
  stationary <- Re(eig$vectors[, i])
  stationary <- stationary / sum(stationary)  # normalize to sum to 1
  # Show stationary distribution
  stationary
}

stat_dist <- t(apply(a_gamma, 3, function(x) getStat(x)))

# How many people have probability less than 5% to be in each state?
sum(stat_dist[, 1] < 0.05)
sum(stat_dist[, 2] < 0.05)
sum(stat_dist[, 3] < 0.05)

# For fixed effects
round(getStat(gamma_group), 2)


# --------------------------------------------------------
# ---------- Figure with 3 Selected Transition Prob ------
# --------------------------------------------------------

# Same three people as in Data Vizualization
u_subj <- unique(emotion_mHMM$subj_id)
n_subj <- length(u_subj)
set.seed(4)
v_sel <- sample(1:n_subj, size=3, replace=F) 
u_subj[v_sel]

# Look at stationary distributions
stat_dist[v_sel, ]

# Select model
m <- 3
model_m <- l_Models[[m]]
gamma_subj <- obtain_gamma(model_m, level = "subject")

# Get random effects

# ----- Plotting -----
sc <- 0.9
pdf("Figures/Transition_3Selected.pdf", width=10*sc, height=4*sc)

# Layout
lmat <- rbind(1:3, 
              4:6)
lo <- layout(lmat, widths = rep(1, 3), heights = c(0.1, 1))
# layout.show(lo)
for(i in 1:3) plotLabel(paste0("Subject ", u_subj[v_sel[i]]), 
                        cex=1.75, xpos=0.55)

# Heatplot/Data
for(i in 1:3) plotHeat(gamma_subj[[v_sel[i]]], cex.axis=1.2, cex.from=1.5, cex.to=1.1)

dev.off()


# --------------------------------------------------------
# ---------- State Decoding ------------------------------
# --------------------------------------------------------

# ---- Load Data -----
emotion_mHMM <- readRDS("Data/data_Rowland2020_wN_DF.RDS")

# ---- Load Model -----
l_Models <- readRDS("Files/RowlandHMMs_rep1.RDS")

# ----- Select People -----
u_subj <- unique(emotion_mHMM$subj_id)
n_subj <- length(u_subj)
set.seed(4)
v_sel <- sample(1:n_subj, size=3, replace=F) 
cols <- brewer.pal(3, "Set1")


# ----- Get Predicted States -----


model_m <- l_Models[[3]]
emotion_states_2st <- vit_mHMM(model_m, emotion_mHMM)
# plot(unique(emotion_states_2st$subj))
head(emotion_states_2st)

u_subj[v_sel]
v_sel

l_decod <- list()
for(i in 1:3) l_decod[[i]] <- emotion_states_2st[emotion_states_2st$subj==u_subj[v_sel[i]], ]

# display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
#                    colorblindFriendly=FALSE)


# ----- Plot Figure -----
cols <- brewer.pal(3, "Set1")
cols2 <- brewer.pal(3, "Dark2") # decoding
cols2 <- c("lightgrey", "darkgrey", "black")

N <- 240 + 40*8 # Number of time points with night-gap
sc <- 0.7

pdf("Figures/Rowland_3subj_1x3panel_withDecoding.pdf", width=12*sc, height=4*sc)

par(mfrow=c(1,3), mar=c(4,4,3.4,0.2))

for(i in 1:3) {
  # Canvas
  plot.new()
  plot.window(xlim=c(0, 600), ylim=c(0,112))
  axis(1)
  title(xlab="Time", line=2.5)
  axis(2, las=2)
  if(i==1) {
    title(ylab="Responses", line=2.5)
  }
  title(main = paste0("Subject ", u_subj[v_sel[i]]), font.main=1, cex=1.2)
  # Decoding
  for(d in 0:559) rect(d, 102, d+1, 112, border=NA, col=cols2[l_decod[[i]]$state[d]])
  # Data
  data_i <- emotion_mHMM[emotion_mHMM$subj_id == u_subj[v_sel[i]], ]
  lines(data_i$happy, col=cols[2])
  lines(data_i$sad, col=cols[1])
  
  # Legend decoding
  if(i==2) {
    mtext("State 1", 3, at=0, adj=0, col=cols2[1], line=-0.25, cex=0.7)
    mtext("State 2", 3, at=235, adj=0, col=cols2[2], line=-0.25, cex=0.7)
    mtext("State 3", 3, at=450, adj=0, col=cols2[3], line=-0.25, cex=0.7)
  }
  
  # Legend
  if(i==2) {
    legend(-60, 98,
           legend=c("Happy"), adj=0,
           text.col=cols[2],
           bty = "n", cex=1.2)
    legend(180, 90,
           legend=c("Sad"),
           text.col=cols[],
           bty = "n", cex=1.2)
  }
} # end i

dev.off()




