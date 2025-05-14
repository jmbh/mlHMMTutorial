## Tutorial: Estimating multilevel Hidden Markov Models

This repository contains the code for a tutorial to estimate multilevel Hidden Markov Models with the R-package mHMMbayes](https://cran.r-project.org/web/packages/mHMMbayes/index.html). In the tutorial, we use the open data from [Rowland & Wenzel (2020)](https://link.springer.com/article/10.1007/s12671-020-01335-4) who shared their materials and data [here on OSF](https://osf.io/jmz2n/). We use the preprocessed version of these data from [Haslbeck et al. (2023)](https://psycnet.apa.org/fulltext/2023-72233-001.html).

The tutorial itself is in the paper:

XXX FILL IN LINK TO PREPRINT HERE XXX

This repository contains the following files and folders:

- `0_Helpers.R` contains helper functions used throughout.
- `1_ProcessData.R` takes the processed data from [Haslbeck et al. (2023)](https://psycnet.apa.org/fulltext/2023-72233-001.html) and adds rows of NAs to correctly represent the "missing" observations during the night. Saves a new dataframe `data_Rowland2020_wN_DF.RDS` in the folder `/Data` which is the dataframe that is being used for all subsequent analyses.
- `2_VisualizeData.R` creates a Figure 8 in the paper, showing the time series of variables Happy and Sad for three randomly selected persons
- `3_EstimateModels.R` Specifies starting values and priors and fits mlHMMs on the data with a sequence of 1:6 states. It runs four chains separately, saves as `RowlandHMMs_repX.RDS` into the folder `/Files`. We run four chains, so that we can assess convergence.
- `4_LabelSwitching.R` Creates trace plots of all parameters, for all chains, for all subjects to assess label switching.
- `5_Convergence.R` uses the four estimated chains to make trace plots (shown in the appendix) and GR-statistics (shown in Table 3) for the models with 1:6 states
- `6_ModelFit.R` takes the fitted model (chain 1), extracts the LL, AIC, and AICc, and compiles them into Table 2 in the paper
- `7_InspectParameters.R` extracts fixed and random effects parameters of emission distributions and transition probability matrices and plots associated figures shown in the paper and the appendix.
- `8_PPCs.R` performs Posterior Predictive Checks (PPCs) on the distribution of within-person means and within-person standard deviations for the models with 1:4 and plots the corresponding figures in the paper and the appendix.
- `9_PseudoResiduals.R` performs a pseudoresidual analysis for the three subjects shown also in Figure 8 for the sequence of models with 1:6 states, and plots associated figures. It also plots the residuals for all persons and all variables, for models with 1:m states, an analysis that is only mentioned and not shown in the paper.




