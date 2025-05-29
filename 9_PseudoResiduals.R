# jonashaslbeck@protonmail.com; May 29th, 2025

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# We perform a Pseudo Residual Analysis for a few variables
# and subjects and plot corresponding figures shown in the
# main text and appendix


# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------

library(mHMMbayes)
library(xtable)
library(plyr)
library(dplyr)
library(RColorBrewer)

library(ggplot2)


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
# ---------- Get Some Info from Data ---------------------
# --------------------------------------------------------

labels <- c("happy", "excited", "relaxed", "satisfied",
            "angry", "anxious", "depressed", "sad")

labels <- colnames(emotion_mHMM)[-1] # Variable names
v_subj_id <- unique(emotion_mHMM$subj_id) # Unique subject IDs



# --------------------------------------------------------
# ---------- Get Numbers for all Variables ---------------
# --------------------------------------------------------

N_subj <- 125
m_R2s <- matrix(NA, N_subj, 8)
m_RMSE <- matrix(NA, N_subj, 8)
m_ARrat <- matrix(NA, N_subj, 8)

l_resid <- list()

for(i in 1:N_subj) {
  
  l_resid_j <- list()
  
  for(j in 1:8) {
    
    # Residual Analysis
    l_resid_j[[j]] <-  res_ji <- GetResid(data = emotion_mHMM,
                                          model = model_m,
                                          j = j,
                                          i = i)
    
    # Save RMSE
    m_RMSE[i, j] <- res_ji$RMSE
    
    ind_NA <- is.na(res_ji$emp)
    
    # Sanity:
    # plot(res_ji$emp, res_ji$model)
    # cor(res_ji$emp, res_ji$model, use="complete.obs")^2
    
    # Compute R2
    m_R2s[i, j] <- 1 - var(res_ji$resid[!ind_NA]) / var(res_ji$emp[!ind_NA])
    
    # AR ratio
    dat_ar_emp <- cbind(res_ji$emp[-1], res_ji$emp[-length(res_ji$emp)])
    ar_emp <- cor(na.omit(ar_ij))[2,1]
    model_res <- res_ji$resid
    model_res[ind_NA] <- NA
    dat_ar_mod <- cbind(model_res[-1], model_res[-length(model_res)])
    ar_res <- cor(na.omit(dat_ar_mod))[2,1]
    m_ARrat[i, j] <- ar_res / ar_emp
  }
  
  l_resid[[i]] <- l_resid_j
  print(i)
}

# ----- Save Results -----
saveRDS(m_R2s, "Files/Residual_R2.RDS")
saveRDS(m_RMSE, "Files/Residual_RMSE.RDS")
saveRDS(m_ARrat, "Files/Residual_ARrat.RDS")

# ---- Values reported in the paper -----
range(round(apply(m_RMSE, 2, median)))
1-range(round(apply(m_ARrat, 2, median), 2))

# ---- Make Figures -----
n_dep <- 8
cols <- brewer.pal(n_dep+1, "Set1")[-6]

df_R2s <- data.frame(m_R2s)
df_ARrat <- data.frame(m_ARrat)
df_RMSE <- data.frame(m_RMSE)
colnames(df_R2s) <- colnames(df_ARrat) <- colnames(df_RMSE)  <- labels

pdf("Figures/ResidualAnalysis_R2s.pdf", width=8, height=5)
boxplot(df_R2s, col=cols, axes=FALSE)
axis(1, labels=labels, at=1:8, las=2)
axis(2, las=2)
grid()
boxplot(df_R2s, col=cols, axes=FALSE, add=TRUE)
title(ylab=expression(R^2))
dev.off()

pdf("Figures/ResidualAnalysis_RMSEs.pdf", width=8, height=5)
boxplot(df_RMSE, col=cols, axes=FALSE)
axis(1, labels=labels, at=1:8, las=2)
axis(2, las=2)
grid()
boxplot(df_RMSE, col=cols, axes=FALSE, add=TRUE)
title(ylab="RMSE")
dev.off()

pdf("Figures/ResidualAnalysis_VARrat.pdf", width=8, height=5)
boxplot(df_ARrat, col=cols, axes=FALSE)
abline(h=1)
axis(1, labels=labels, at=1:8, las=2)
axis(2, las=2)
grid()
boxplot(df_ARrat, col=cols, axes=FALSE, add=TRUE)
title(ylab="AR(residual)  / AR(data)")
dev.off()

# ----- Some sanity checking -----
which(df_R2s[, 6] < -0.8)
v_subj_id[87]

res_ji <- GetResid(data = emotion_mHMM,
                   model = model_m,
                   j = 6,
                   i = 87)

data_s <- emotion_mHMM[emotion_mHMM$subj_id == v_subj_id[87], ]
plot(res_ji$emp, type="l", ylim=c(0,8))
lines(res_ji$model, col="red")
legend("topright", legend=c("data", "HMM prediction"), text.col=c("black", "red"), bty="n")


# --------------------------------------------------------
# ---------- Figure: 3 Persons, m States -----------------
# --------------------------------------------------------

# For m=3: Figure 10 in the paper
# m=2,4 are shown in the appendix
# m=1,5,6 are also plotted here, but not shown in the paper

# Get same 3 random subjects as in the data inspection figure
set.seed(4)
v_sel <- sample(1:n_subj, size=3, replace=F) # Sample 4 people at random

# Loop through m
for(m in 1:6) {
  
  sc <- 1.05
  pdf(paste0("Figures/Fig_ResidualsPaper_3subj_m=", m, ".pdf"), width=7.5*sc, height=6.5*sc)
  
  # Select Model
  model_m <- l_Models[[m]]
  
  # Layout
  lmat <- matrix(6:14, 3, 3, byrow = TRUE)
  lmat <- rbind(c(1:2, 0), lmat)
  lmat <- cbind(c(0, 3:5), lmat)
  lo <- layout(mat = lmat, widths = c(0.15, 1, 1, .3), heights = c(0.15, rep(1, 4)))
  # layout.show(lo)
  
  # Plot Labels
  plotLabel("Data + Predictions", cex=1.6)
  plotLabel("Pseudo Residuals", cex=1.6)
  for(s in 1:3) plotLabel(paste0("Subject ", v_subj_id[v_sel[s]]), srt=90)
  
  for(s in 1:3) {
    
    # Get predictions, residuals
    res_1i <- GetResid(data = emotion_mHMM, 
                       model = model_m,
                       j = 1, # Variable
                       i = v_sel[s]) # This function is in 0_Helpers.R
    
    # Plot data + Predictions
    par(mar=c(4,2,1,1))
    plot.new()
    plot.window(xlim=c(0, 600), ylim=c(0,100))
    axis(1)
    axis(2, las=2)
    lines(res_1i$emp, lwd=1.5)
    lines(res_1i$model, col="orange", lty=2, lwd=1.5)
    if(s==1) legend("bottomright", legend=c("Data 'Happy'", "Prediction"),
                    text.col=c("black", "orange"), lty=1:2, col=c("black", "orange"),
                    bty="n")
    
    # Plot residuals
    PlotRes(res_1i, 
            j = 1, 
            i = v_sel[s], 
            layout = FALSE, 
            title = FALSE) # This function is in 0_Helpers.R
    
  } # end for: subj
  
  dev.off()
  
  
} # end for:m


# --------------------------------------------------------
# ---------- All Persons, m States -----------------------
# --------------------------------------------------------

# This is only mentioned, not shown in the paper, only included in this repo

# Here we run this over states 1:4, which are the ones we 
# consider after the convergence check

for(m in 1:4) {
  
  # Select model with m states
  model_m <- l_Models[[m]]
  
  sc <- 1.3
  pdf(paste0("Figures/ResidualAnalysis_all_m=", m, ".pdf"), width=6*sc, height=7*sc)
  
  N_subj <- length(v_subj_id)
  
  # Storage for RMSEA across people and items
  m_RMSE <- matrix(NA, N_subj, 8)
  
  # Make layout
  lmat <- matrix(1:16, ncol=4, byrow = TRUE)
  lo <- layout(mat=lmat, widths = rep(c(1, .5), 2))
  
  for(i in 1:N_subj) {
    for(j in 1:8) {
      
      res_ji <- GetResid(data = emotion_mHMM,
                         model = model_m,
                         j = j,
                         i = i)
      m_RMSE[i, j] <- res_ji$RMSE
      
      PlotRes(res_ji, j=j, i=i, layout=FALSE)
      
    }
    print(i)
  }
  
  dev.off()
  
} # end for: m



# --------------------------------------------------------
# ---------- Look at Residual Correlations ---------------
# --------------------------------------------------------


# ----- Get predicted states for person -----
# Needed to condition on estimated states
l_Models <- readRDS("Files/RowlandHMMs_rep1.RDS")
model_m <- l_Models[[3]]
emotion_states_2st <- vit_mHMM(model_m, emotion_mHMM)

# ---------- Get residual matrix & cors for every person --------

l_ResDat <- list()
l_ResCorm <- list()
l_ResCorm1 <- l_ResCorm2 <- l_ResCorm3 <- list()

for(i in 1:N_subj) {
  m_resid <- matrix(NA, length(l_resid[[i]][[1]]$resid), 8)
  for(j in 1:8) m_resid[, j] <- l_resid[[i]][[j]]$resid
  l_ResDat[[i]] <- m_data
  l_ResCorm[[i]] <- cor(m_resid, use="complete.obs")
  
  # Compute Residual covs separately for states
  states_i <- emotion_states_2st$state[emotion_states_2st$subj==v_subj_id[i]] 
  if(sum(states_i==1) > 1) l_ResCorm1[[i]] <- cor( m_resid[states_i==1, ], use="complete.obs") else l_ResCorm1[[i]] <- matrix(NA, 8, 8)
  if(sum(states_i==2) > 1) l_ResCorm2[[i]] <- cor( m_resid[states_i==2, ], use="complete.obs") else l_ResCorm2[[i]] <- matrix(NA, 8, 8)
  if(sum(states_i==3) > 1) l_ResCorm3[[i]] <- cor( m_resid[states_i==3, ], use="complete.obs") else l_ResCorm3[[i]] <- matrix(NA, 8, 8)
}

a_ResCorm <- simplify2array(l_ResCorm)

a_ResCorm1 <- simplify2array(l_ResCorm1)
a_ResCorm2 <- simplify2array(l_ResCorm2)
a_ResCorm3 <- simplify2array(l_ResCorm3)



# ---------- Overall: Mean & Inspect --------

# Compute
ResCorm_mean <- apply(a_ResCorm, 1:2, mean)
round(ResCorm_mean, 2)

plotHeat <- function(corm, title=NULL) {
  
  # Visualize
  cor_mat <- corm
  rownames(cor_mat) <- labels
  colnames(cor_mat) <- labels
  
  # Convert to long format
  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Var1", "Var2", "Correlation")
  
  # 2. Filter to upper triangle only (excluding diagonal)
  cor_df <- cor_df[as.numeric(factor(cor_df$Var1, levels = labels)) <
                     as.numeric(factor(cor_df$Var2, levels = labels)), ]
  
  # 3. Add text labels
  cor_df$label <- sprintf("%.2f", cor_df$Correlation)
  
  # Plot
  ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile() +
    geom_text(aes(label = label), size = 4) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1, 1)) +
    theme_minimal() +
    labs(x = NULL, y = NULL, title=title) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}
p1 <- plotHeat(ResCorm_mean, title="All States")


# Alright, quite some residuals left

# SDs
ResCorm_SD <- apply(a_ResCorm, 1:2, sd)
round(ResCorm_SD, 2) # huge between-person variation in residual correlations


# ---------- State-Specific: Mean & Inspect --------

ResCorm1_mean <- apply(a_ResCorm1, 1:2, mean, na.rm=T)
round(ResCorm1_mean, 2)
p2 <- plotHeat(ResCorm1_mean, title="State = 1")

ResCorm2_mean <- apply(a_ResCorm2, 1:2, mean, na.rm=T)
round(ResCorm2_mean, 2)
p3 <- plotHeat(ResCorm2_mean, title="State = 1")

ResCorm3_mean <- apply(a_ResCorm3, 1:2, mean, na.rm=T)
round(ResCorm3_mean, 2)
p3 <- plotHeat(ResCorm3_mean, title="State = 1")

# Plot layout
library(patchwork)
(p1 | p2) / (p3 | p4)
















