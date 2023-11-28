Here is a supplementary program and dataset for the manuscript "Flexible Modeling of Nonstationary Extremal Dependence Using Spatially-Fused LASSO and Ridge Penalties".


R files:

1. (*)Modeling.R: Modeling extremal dependence with given data;
2. (*)Summary.R: Summary of the fitted model with related plots and tables;
3. Simulate_data.R: A simple simulation of max-stable processes with nonstationary extremal dependence;
4. Algorithm1.R: Function for the Algorithm 1 described in the paper;
5. Merge_subr.R: Function for the subregion merging process;
6. Lambda_tuning.R: Function for the lambda-tuning;
7. Fit.R: Some fitting function for convenience, based on the optim;
8. Objectives.R: Pairwise likelihood function for the max-stable processes (and also for the inverted MSPs);
9. Utils.R: Some auxiliary functions.
(*): main scripts


Rdata files:

1. NepalExtended.Rdata: The gridded monthly maximum temperature dataset from Nepal and its surrounding Himalayan and sub-Himalayan regions used in data application. The file includes the parameter estimation results of margin structure estimation (GEV parameters) using the Max-and-Smooth method: postman.mu (location), postman.sigma (scale), postman.xi (shape).

2. Simulated.Rdata: A simulated data with Simulate_data.R, including the coordinate and true (dependence) parameter information. The true partition is the partition P1 mentioned in the paper.

Remarks:

1. The current program only works for grid data, but an extension to general lattice data is available.
2. Nonstationarity is assumed for the input data.
3. Input data must be renormalized in unit Fréchet scale for the max-stable processes.
4. Some difficulties in range estimation may emerge.

How to Start:

1. Input the dataset (simulation or extended Nepal data) by changing the variable "dat". Hyperparameters are given according to the chosen dataset.
2. Compile Modeling.R for modeling and Summary.R for assessment with plots and tables, which are related to Figures 4-6 and Table 3 in the paper.
