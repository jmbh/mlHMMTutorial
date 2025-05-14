# jonashaslbeck@protonmail.com; Feb 5th, 2025

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Helper functions


# --------------------------------------------------------
# ---------- Plotting Labels in into Canvas --------------
# --------------------------------------------------------

plotLabel <- function(text, cex=1.4, srt=0, xpos=0.5) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(xpos, 0.5, text, cex=cex, srt=srt)
}


# --------------------------------------------------------
# ---------- Compute Gelman-Ruben Statistic --------------
# --------------------------------------------------------

f_GR <- function(model_list, digits=2, burnin=500) {
  
  # Get m
  nChain <- length(model_list) # number of Chains
  m <- model_list[[1]]$input$m
  
  # Transition Matrix
  iter <- m*(m-1)
  v_RG_gamma <- rep(NA, m*(m-1))
  if(m>1) {
    for(j in 1:iter) {
      l_obj <- lapply(model_list, function(x) mcmc(x$gamma_int_bar[-(1:burnin), j]))
      mcmc_chains <- mcmc.list(l_obj) #make this a list input (must be possible)
      gelman_rubin_stat <- gelman.diag(mcmc_chains)
      v_RG_gamma[j] <- gelman_rubin_stat$psrf[1,1]
    }
  }
  m_RG_gamma <- matrix(v_RG_gamma, m, m-1, byrow = TRUE)
  
  # Emission Distribution
  RG_emiss <- matrix(NA, n_dep, m)
  for(p in 1:n_dep) {
    for(ms in 1:m) {
      l_obj <- lapply(model_list, function(x) mcmc(x$emiss_mu_bar[[p]][-(1:burnin), ms]))
      mcmc_chains <- mcmc.list(l_obj)
      gelman_rubin_stat <- gelman.diag(mcmc_chains)
      RG_emiss[p, ms] <- gelman_rubin_stat$psrf[1,1]
    }
  }
  
  # Return
  outlist <- list("m_RG_gamma" = round(m_RG_gamma, digits=digits),
                  "RG_emiss" = round(RG_emiss, digits=digits))
  return(outlist)
  
} # eoF


# --------------------------------------------------------
# ---------- Functions to Compute Pseudo-Residuals -------
# --------------------------------------------------------

GetResid <- function(data, # empirical data
                     model, # fitted mHMM model object
                     j, # Desired variable
                     i # Desired subject
) {
  
  # Aux
  v_subj_id <- unique(data$subj_id)
  
  # --- Emission Distributions ---
  emiss_subject <- obtain_emiss(model, level = "subject") # For all
  mu_ji <- emiss_subject[[labels[j]]][[i]][, 1] # For specified variable & subject; col=1 for means
  
  # Get most likely state sequence for each subject
  if(model$input$m > 1) {
    model_m_stateSeq <- suppressMessages(vit_mHMM(model, emotion_mHMM))
    State_i <- model_m_stateSeq$state[model_m_stateSeq$subj==v_subj_id[i]] # state seq for person i
  } else {
    State_i <- rep(1, 240+40*8) # if m=1, we know that state is 1 always
  }
  
  # Compute predictions for fixed j & i
  Pred_ji <- as.numeric(mu_ji[State_i])
  
  # Get Empirical data
  X_ji <- data[[labels[j]]][data$subj_id==v_subj_id[i]]
  
  # Compute residual
  res_ji <- X_ji - Pred_ji
  
  # Compute MSE
  RMSE_ji <- sqrt(mean(na.omit(res_ji)^2))
  
  # Return
  outlist <- list("emp" = X_ji,
                  "model" = Pred_ji,
                  "resid" = res_ji,
                  "RMSE" = RMSE_ji)
  return(outlist)
  
} # eoF


# --------------------------------------------------------
# ---------- Functions Visualize Pseudo-Residuals --------
# --------------------------------------------------------

PlotRes <- function(res_ji, j, i, layout=TRUE, title=TRUE) {
  
  # Compute AR
  ar_ji <- round(cor(res_ji$resid[-1], res_ji$resid[-(240+8*40)], use="complete.obs"), 3)
  # Fit linear trend
  time <- 1:(240+8*40)
  trend_model <- lm(res_ji$resid ~ time)
  trend_sum <- summary(trend_model)
  trend_coef <- round(trend_sum$coefficients[2, 1], 3)
  trend_pval <- round(trend_sum$coefficients[2, 4], 3)
  
  # Layout
  if(layout) layout(matrix(1:2, ncol=2), widths = c(1, .5))
  # LinePlot
  par(mar=c(4,3,2,1))
  plot.new()
  plot.window(xlim=c(1, 560), ylim=c(-110, 110))
  axis(1)
  axis(2, las=2)
  abline(h=0, lty=1, col="lightgrey", lwd=2)
  points(res_ji$resid, pch=20)
  abline(trend_model, lwd=2, col="steelblue")
  text(5,-70, paste0("Lin Trend: Slope = ", trend_coef, "; pval = ", trend_pval), col="steelblue", adj=0)
  text(5,-90, paste0("Lag-1 AR = ", ar_ji), col="tomato", adj=0)
  text(4, 80, paste0("RMSE = ", round(res_ji$RMSE, 2)), adj=0)
  if(title) title(main=paste0(labels[j], " [subj = ", v_subj_id[i], "]"), font.main=1)
  # Marginal
  par(mar=c(4,0,2,1))
  hist_data <- hist(res_ji$resid, plot = FALSE, breaks=seq(-110, 110, length=20))
  barplot(hist_data$counts,
          horiz = TRUE,  # Horizontal bars
          names.arg = NULL,
          axes=FALSE)
  x_seq <- seq(-100, 100, length=1000)
  gauss_den <- dnorm(x_seq,
                     mean = mean(res_ji$resid, na.rm = TRUE),
                     sd = sd(res_ji$resid, na.rm = TRUE))
  scaled_den <- (gauss_den/max(gauss_den) ) * max(hist_data$counts)
  lines(scaled_den, seq(0, 24, length=1000), col="grey") # Not entiresure where those 22 come from
  
} # eoF

# --------------------------------------------------------
# ---------- Extract Transition Estimates & Plot ---------
# --------------------------------------------------------


PlotTrans <- function(object = l_Models,
                      m = m,
                      N_disp = 50) # number of random effects shown in the figure
  
{
  
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
  # Layout
  lo <- layout(matrix(1:2, ncol=2), widths = c(.4, .6))
  # layout.show(lo)
  
  ### Left Panel: Transition Matrix
  plt_h <- plotHeat(gamma_group, cex.axis=1.2, cex.from=1.5, cex.to=1.5)
  
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
  cnt <- 1
  for(m1 in 1:m) {
    for(m2 in 1:m) {
      points(rep(bp[cnt], N_disp), a_gamma[m1, m2, 1:N_disp],
             pch=21, col="black",
             bg=as.character(plt_h[m:1,])[cnt], cex=0.75)
      cnt <- cnt + 1
    }
  }
  
  
  
} # eoF



# --------------------------------------------------------
# ---------- Extract Emission Estimates & Plot -----------
# --------------------------------------------------------

PlotEmiss <- function(object, m, parameter="Means") {
  
  # Basic info (hardcoded)
  n_dep <- 8 
  N <- 125
  
  model_m <- object[[m]]
  
  # Parameter indicator
  indpar <- ifelse(parameter=="Means", 1, 2)
  
  # ------ Getting the Parameters ------
  emiss_group <- obtain_emiss(model_m)
  labels <- names(emiss_group)
  
  ## Fixed effects
  m_fixed <- matrix(NA, n_dep, m)
  for(i in 1:n_dep) m_fixed[i, ] <- emiss_group[[i]][, indpar] # 1=mean; 2=SD; later make variable
  colnames(m_fixed) <- paste("K = ", 1:m)
  rownames(m_fixed) <- labels
  # Order effects
  ord <- order(m_fixed[1,], decreasing = TRUE)
  m_fixed_ord <- m_fixed[, ord]
  
  ## Random Effects
  emiss_subject <- obtain_emiss(model_m, level = "subject")
  a_REs <- array(NA, dim=c(n_dep, m, N))
  for(j in 1:N) for(i in 1:n_dep) a_REs[i, , j] <- emiss_subject[[i]][[j]][, indpar] # 1=mean; 2=SD; later make variable
  # Order REs
  a_REs_ord <- a_REs[, ord, ]
  
  # ------ Plotting ------
  
  ## Plotting basics
  library(RColorBrewer)
  cols <- brewer.pal(n_dep+1, "Set1")[-6]
  
  pdf(paste0("Figures/Emission_", parameter, "_m=", m, ".pdf"), width = 7, height = 4)
  
  par(mar=c(3,3,2,1))
  old_par <- par(no.readonly = TRUE)
  par(mgp = c(3, 0.25, 0)) # moving x-axis labels up a bit
  
  ## Get started
  bp <- barplot(t(m_fixed_ord),
                beside=TRUE,
                space = c(0.2, 1),
                ylim=c(0, 100),
                col = rep(cols, each=m),
                names.arg = rep(1:m, n_dep),
                cex.names=0.7, axes=FALSE)
  # Emotion labels
  em_m <- colMeans(bp)
  for(i in 1:n_dep) text(em_m[i], 98, labels[i], adj=0.5, col=cols[i], cex=0.85)
  # Restore old par
  # Axes/labels/etc
  title(xlab="States", line=1.5)
  
  ## Add random effects
  # REs are available only for means
  if(indpar==1) {
    
    N_disp <- 50 #N
    for(i in 1:n_dep) {
      # Plot lines between dots
      for(s in 1:(m-1)) {
        segments(x0 = bp[s, i],
                 y0 = a_REs_ord[i, s, 1:N_disp],
                 x1 = bp[s+1, i],
                 y1 = a_REs_ord[i, s+1, 1:N_disp], col=alpha("grey", alpha=0.75), lwd=0.5)
      }
      
      # Plot those dots
      for(s in 1:m) points(rep(bp[s, i], N_disp), a_REs_ord[i, s, 1:N_disp],
                           col="black", bg=cols[i], pch=21, cex=0.5)
      
    } # end for: vars
  } # end if: means?
  
  par(old_par)
  axis(2, las=2, seq(0, 100, length=6), at=seq(0, 1, length=6))
  
  dev.off()
  
  
} # eoF


# --------------------------------------------------------
# ---------- Heatplot [for Transition Matrices] ----------
# --------------------------------------------------------

plotHeat <- function(phi,
                     k,
                     main="",
                     labels=NULL,
                     las.x=1,
                     cex.axis=1.5,
                     cex.from = 1.5,
                     cex.to=1.5,
                     cex.val=1) {
  
  
  # Save the current par settings
  old_par <- par(no.readonly = TRUE)
  
  # -- Aux Variables --
  p <- ncol(phi)
  
  # -- Make color gradient --
  color.gradient <- function(x, colors=c("#E41A1C", "white", "#377EB8"), colsteps=201) {
    return( grDevices::colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  x <- 1:201
  grad <- color.gradient(x)
  
  # Make canvas
  par(mar=c(2,2.5,2.5,1)*2)
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0, 1))
  
  # Auxiliary plotting variables
  sfm <- 1/(p*2)
  seq_mp_x <- seq(0, 1, length=p+1)[-(p+1)] + sfm
  
  # Plot Axes & Axis labels
  # xy_labels <- paste0("S ", 1:p)
  # Adjust mgp for custom axis labels
  par(mgp = c(3, 0.35, 0))  # Move tick labels closer to the plot
  
  y_labels <- sapply(p:1, function(i) as.expression(bquote(S[.(i)])))
  x_labels <- sapply(1:p, function(i) as.expression(bquote(S[.(i)])))
  
  axis(3, labels = x_labels, at=seq_mp_x, cex.axis=cex.axis, tick=FALSE)
  axis(2, labels = y_labels, at=seq_mp_x, las=2, cex.axis=cex.axis, tick=FALSE)
  title(main, font.main=1)
  
  title(ylab="From", cex.lab=cex.from, line=2)
  mtext("To", side=3, cex=cex.to, line=2)
  
  phi_col <- matrix(NA, p, p)
  
  # Plot Data
  for(i in 1:p) {
    for(j in 1:p) {
      
      # Get color
      phi_ij <- phi[p:1, ][j, i]
      if(phi_ij < -1) {
        col_ij <- grad[1]
      } else if(phi_ij > 1 ) {
        col_ij <- grad[201]
      } else {
        col_ij <- grad[phi[p:1, ][j, i] * 100 + 101]
        phi_col[j,i] <- col_ij
      }
      
      # Plot box
      rect(xleft = seq_mp_x[i]-sfm,
           ybottom = seq_mp_x[j]-sfm,
           xright = seq_mp_x[i]+sfm,
           ytop = seq_mp_x[j]+sfm,
           col = col_ij)
      # Plot text
      text(seq_mp_x[i], seq_mp_x[j], round(phi_ij , 2), cex=cex.val, col="black")
    }
  }
  
  # Reset par settings back to original
  par(mgp = c(3, 1, 0))
  
  
  # Return colors
  return(phi_col)
  
} # eoF







