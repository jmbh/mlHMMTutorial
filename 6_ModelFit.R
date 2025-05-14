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


# --------------------------------------------------------
# ---------- Load Model Object ---------------------------
# --------------------------------------------------------

l_Models <- readRDS("Files/RowlandHMMs_rep1.RDS")


# --------------------------------------------------------
# ---------- Get LL/AIC/AICc -----------------------------
# --------------------------------------------------------

# ---------- Look at Adj AIC ---------

l_Models[[1]]
l_Models[[2]]
l_Models[[3]]
l_Models[[4]]
l_Models[[5]]
l_Models[[6]]

v_LL <- c(-6249.804, -5837.17, -5680.119, -5574.184 , -5514.328, -5482.26)
v_AIC <- c(12531.61, 11742.34, 11468.24, 11300.37, 11228.66, 11216.52)
v_AICc <- c(12532.61, 11746.87, 11480, 11324.6, 11272.67, 11290.43) 

# Negative LL
v_nLL <- - v_LL
# Scale to max IACc
v_nLL_sc <- (v_nLL/max(v_nLL)) * max(v_AICc)

# ----- Make Table -----
tab_fit <- data.frame(round(rbind(v_nLL, v_AIC, v_AICc)))
colnames(tab_fit) <- paste0("M = ", 1:6)
rownames(tab_fit) <- c("-LL", "AIC", "AICc")

# Drop M = 5,6 models, because we excluded them due to convergence issues
tab_fit <- tab_fit[, 1:4]

xtbl <- xtable(tab_fit, digits=rep(0, 4+1))
xtbl # This table is shown in the paper








