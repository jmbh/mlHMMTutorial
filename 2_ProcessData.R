# jonashaslbeck@protonmail.com; Feb 5th, 2025

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# We add rows of NAs to the data to represent the night gap
# We do this both to show the night gap in the data visualization
# and because this is how we give the mlHMM estimation function
# the information that this is how many data points are "missing"
# over night

# --------------------------------------------------------
# ---------- Load Data -----------------------------------
# --------------------------------------------------------

# This is the data as shared by Haslbeck et al. (2023) in Emotion
data <- readRDS("Data/data_Rowland2020.RDS")


# --------------------------------------------------------
# ---------- Add NAs for Night Gap -----------------------
# --------------------------------------------------------

# ----- Determine Number Missing Datapoints in Night -----

# There were six measurements between 10h and 20h
# There was some randomization, but measurements were at least 45min apart
# This means that the average duration between measurements (i.e., the lag) is:
lag = (20-10)/6
lag
# The unmeasured "night time" is
night <- 4+10
night
# And so we need to have
night / lag
# missing values during the night; we round this down to 8


# ----- Insert them in the Data -----
# Unique subjects
u_subjid <- sort(unique(data$subj_id), decreasing = FALSE)
n_subj <- length(u_subjid)

# Storage
l_data_wN <- list()

# Loop over subjects
for(i in 1:n_subj) {
  # Subset
  data_i <- as.matrix(data[data$subj_id==u_subjid[i], ])
  
  # Define Night block
  night_block <- matrix(NA, 8, 12)
  night_block[, 1] <- u_subjid[i]
  
  # Split by day (easy, because we have the full rows for all individuals)
  l_data_i_wN <- list()
  for(day in 1:40) {
    l_data_i_wN[[day]] <- rbind(data_i[data_i[, 2]==day, ],
                                night_block)
  }
  # Combine
  data_i_wN <- do.call(rbind, l_data_i_wN)
  l_data_wN[[i]] <- data_i_wN
} # end loop

# Combine dataset
data_wN <- as.data.frame(do.call(rbind, l_data_wN))
colnames(data_wN) <- colnames(data)

head(data_wN)

# ----- Rearrange dataframe -----
emotion_mHMM <- data.frame(subj_id = data_wN$subj_id,
                           happy = data_wN$happy,
                           excited = data_wN$excited,
                           relaxed = data_wN$relaxed,
                           satisfied = data_wN$satisfied,
                           angry = data_wN$angry,
                           anxious = data_wN$anxious,
                           depressed = data_wN$depressed,
                           sad = data_wN$sad)

# --------------------------------------------------------
# ---------- Save ----------------------------------------
# --------------------------------------------------------

saveRDS(emotion_mHMM, "Data/data_Rowland2020_wN_DF.RDS")

# These data are used for all subsequent analyses



