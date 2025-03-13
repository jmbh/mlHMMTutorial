# jonashaslbeck@protonmail.com; March 13th, 2025

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# We check whether label swtiching occurs for, first, one
# person
# To do this, we need to check for emission mean sequences whether they are switching


# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------

# install.packages("RColorBrewer")

library(devtools)
library(mHMMbayes)

library(scales)
library(RColorBrewer)


# --------------------------------------------------------
# ---------- Load Estimated Models -----------------------
# --------------------------------------------------------

l_Models <- readRDS("Files/RowlandHMMs_wN_Oct_1_24_m1to6rep1.RDS")


# --------------------------------------------------------
# ---------- Get Posterior Smaples of Person 1 -----------
# --------------------------------------------------------

m <- 3
model_m <- l_Models[[m]]

# --- Gamma ---

# Random effects
dim(model_m$gamma_int_subj[[1]])

# --- Emissions ---
# Fixed effects
dim(model_m$emiss_mu_bar[[1]])
dim(model_m$emiss_mu_bar[[1]])
head(model_m$emiss_mu_bar[[1]])
model_m$emiss_mu_bar$happy

# Random effects; subject 1
dim(model_m$PD_subj[[1]]$cont_emiss)
head(model_m$PD_subj[[1]]$cont_emiss)
model_m$PD_subj[[1]]$cont_emiss[, 1]



# --------------------------------------------------------
# ---------- Trace Plots for all people [all variables] ------------------
# --------------------------------------------------------

# Labels
labels <- c("happy", "excited", "relaxed", "satisfied",
            "angry", "anxious", "depressed", "sad")

# Pick colors
cols <- alpha(brewer.pal(4, "Set1"), alpha=1)

# Loop over M
for(m in 2:4) {

  model_m <- l_Models[[m]]

  # ----- Plotting -----
  pdf(paste0("Figures/Fig_Tutorial_CheckLabelSwitching_ALL_vars_ALL_pers_M=", m, ".pdf"), width=24, height=9)

  cnt <- 1
  for(page in 1:25) {

    # Create layout
    lmat <- rbind(c(0, 1:8),
                  c(9, 14:21),
                  c(10, 22:29),
                  c(11, 30:37),
                  c(12, 38:45),
                  c(13, 46:53))
    lo <- layout(mat=lmat, widths = c(0.2, rep(1, 8)), heights = c(0.2, rep(1, 5)))
    # layout.show(lo)

    # Plot labels
    for(i in 1:8) plotLabel(labels[i])
    for(i in 1:5) plotLabel(paste0("Person ", cnt+i-1), srt=90)


    par(mar=c(2,2,1,1))

    # Five persons per page
    for(i in 1:5) {

      # Loop over 8 dependent variables
      for(v in 1:8) {

        # Canvas
        plot.new()
        plot.window(xlim=c(0, 2000), ylim=c(-10, 110))
        axis(1)
        axis(2, las=2)
        grid()
        if(i==1) legend("bottom", legend=paste0("State ", 1:m), text.col=cols[1:m], bty="n", horiz = TRUE)

        # Data
        for(j in 1:m) lines(model_m$PD_subj[[cnt]]$cont_emiss[, v*m-m +j], type="l", col=cols[j], lwd=1)

      } # end for: variables

      # Update counter
      cnt <- cnt + 1
    }

  } # end for: page

  dev.off()

} # end for: k



# --------------------------------------------------------
# ---------- Make Figure: Sad/Happy k=3, for 5 people ----
# --------------------------------------------------------


# ----- Plotting -----
pdf("Figures/Fig_Tutorial_CheckLabelSwitching_5pers.pdf", width=6, height=9)

# Create layout
lmat <- rbind(c(0, 1, 2),
              c(3, 8, 9),
              c(4, 10, 11),
              c(5, 12, 13),
              c(6, 14, 15),
              c(7, 16, 17))
lo <- layout(mat=lmat, widths = c(0.2, 1, 1), heights = c(0.2, rep(1, 5)))
# layout.show(lo)

# Plot labels
plotLabel("Mean Happy")
plotLabel("Mean Sad")
for(i in 1:5) plotLabel(paste0("Person ", cnt+i-1), srt=90)


par(mar=c(2,2,1,1))
for(i in 1:5) {

  plot.new()
  plot.window(xlim=c(0, 2000), ylim=c(-10, 110))
  axis(1)
  axis(2, las=2)
  grid()
  if(i==1) legend("bottom", legend=paste0("State ", 1:3), text.col=cols, bty="n", horiz = TRUE)
  # Happy
  for(j in 1:3) lines(model_m$PD_subj[[i]]$cont_emiss[, j], type="l", col=cols[j], lwd=1)

  plot.new()
  plot.window(xlim=c(1000, 2000), ylim=c(-10, 110))
  axis(1)
  axis(2, las=2)
  grid()
  # Sad
  for(j in 22:24) lines(1000:2000, model_m$PD_subj[[i]]$cont_emiss[1000:2000, j], type="l", col=cols[j-21], lwd=1)

}

dev.off()


# --------------------------------------------------------
# ---------- Same Figure for all people ------------------
# --------------------------------------------------------

# ----- Plotting -----
pdf("Figures/Fig_Tutorial_CheckLabelSwitching_ALL_pers.pdf", width=6, height=9)

cnt <- 1
for(page in 1:25) {

  # Create layout
  lmat <- rbind(c(0, 1, 2),
                c(3, 8, 9),
                c(4, 10, 11),
                c(5, 12, 13),
                c(6, 14, 15),
                c(7, 16, 17))
  lo <- layout(mat=lmat, widths = c(0.2, 1, 1), heights = c(0.2, rep(1, 5)))
  # layout.show(lo)

  # Plot labels
  plotLabel("Mean Happy")
  plotLabel("Mean Sad")
  for(i in 1:5) plotLabel(paste0("Person ", cnt+i-1), srt=90)

  par(mar=c(2,2,1,1))

  for(i in 1:5) {

    plot.new()
    plot.window(xlim=c(0, 2000), ylim=c(-10, 110))
    axis(1)
    axis(2, las=2)
    grid()
    if(i==1) legend("bottom", legend=paste0("State ", 1:3), text.col=cols, bty="n", horiz = TRUE)
    # Happy
    for(j in 1:3) lines(model_m$PD_subj[[cnt]]$cont_emiss[, j], type="l", col=cols[j], lwd=1)

    plot.new()
    plot.window(xlim=c(1000, 2000), ylim=c(-10, 110))
    axis(1)
    axis(2, las=2)
    grid()
    # Sad
    for(j in 22:24) lines(1000:2000, model_m$PD_subj[[cnt]]$cont_emiss[1000:2000, j], type="l", col=cols[j-21], lwd=1)

    # Update counter
    cnt <- cnt + 1
  }

} # end for: page

dev.off()














