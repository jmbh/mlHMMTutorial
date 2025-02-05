# jonashaslbeck@protonmail.com; Feb 5th, 2025

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# We visualize the time series of the variables "Happy"
# and "Sad" for three different subjects


# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------

library(RColorBrewer)


# --------------------------------------------------------
# ---------- Load Data -----------------------------------
# --------------------------------------------------------

emotion_mHMM <- readRDS("Data/data_Rowland2020_wN_DF.RDS")


# --------------------------------------------------------
# ---------- Visualize Data: 1x3 -------------------------
# --------------------------------------------------------

# Sample 3 people at random
u_subj <- unique(emotion_mHMM$subj_id)
n_subj <- length(u_subj)
set.seed(4)
v_sel <- sample(1:n_subj, size=3, replace=F) 
cols <- brewer.pal(3, "Set1")

N <- 240 + 40*8 # Number of time points with night-gap

sc <- 0.7
pdf("Figures/Rowland_3subj_1x3panel.pdf", width=12*sc, height=4*sc)

par(mfrow=c(1,3), mar=c(4,4,3.4,0.2))

for(i in 1:3) {
  # Canvas
  plot.new()
  plot.window(xlim=c(0, 600), ylim=c(0,100))
  axis(1)
  title(xlab="Time", line=2.5)
  axis(2, las=2)
  if(i==1) {
    title(ylab="Responses", line=2.5)
  }
  title(main = paste0("Subject ", u_subj[v_sel[i]]), font.main=1, cex=1.2)
  # Data
  data_i <- emotion_mHMM[emotion_mHMM$subj_id == u_subj[v_sel[i]], ]
  
  
  lines(data_i$happy, col=cols[2])
  lines(data_i$sad, col=cols[1])
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

