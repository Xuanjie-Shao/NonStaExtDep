Here is supplementary program and dataset for the manuscript "Flexible Modeling of Nonstationary Extremal Dependence Using Spatially-Fused LASSO and Ridge Penalties".


R files:
1. (*)Modeling.R: modeling extremal dependence with given data;
2. (*)Summary.R: summary of the fitted model with related plots and tables;
3. Simulate_data.R: a simple simulation of max-stable processes with nonstationary extremal dependence;
4. Algorithm1.R: Function for the Algorithm 1 described in the paper;
5. Merge_subr.R: Function for the subregion merging process;
6. Lambda_tuning.R: Function for the lambda-tuning;
7. Fit.R: some fitting function for convenience, based on the optim;
8. Objectives.R: pairwise likelihood function for the max-stable processes (and also for the inverted MSPs);
9. Utils.R: some auxiliary functions.
(*): main scripts


Rdata files:
1. NepalExtended.Rdata: extended Nepal data used in data application, which includes the parameter estimation results of margin structure estimation (GEV parameters) using Max-and-Smooth method: postman.mu (location), postman.sigma (scale), postman.xi (shape).
2. Simulated.Rdata: a simulated data with Simulate_data.R, including the coordinate and true (dependence) parameter information. The true partition is the partition P1 mentioned in the paper.


Remarks:
1. Current program only works for grid data, but extension to general lattice data is available.
2. Nonstationarity is assumed for the input data.
3. Input data must be renormalised in unit Fréchet scale for the max-stable processes.
4. Some difficulties in range estimation may emerge.


How to Start:
1. Provide the work directory before compiling each R file.
2. Input the dataset (simulation or extended Nepal data) by changing the variable "dat". Hyperparameters are given according to the chosen dataset.
3. Compile Modeling.R for modeling and Summary.R for assessment with plots and tables, which are related to the Figures 4-6 and Table 3 in the paper.



