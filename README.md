# Flexible Modeling of Nonstationary Extremal Dependence Using Spatially-Fused LASSO and Ridge Penalties
Functions for modeling nonstationary extremal dependence using locally-statioanry max-stable processes with LASSO and ridge penalization.
The provided code is in support of the paper:
<ul> 
          <li> Shao, X., Hazra, A., Richards, J., and Huser, R. (2023+). Flexible modeling of non-stationary extremal dependence using spatially-fused LASSO and ridge penalties. <u><a href="https://arxiv.org/abs/2210.05792" download>ArXiv</a></u> </li>
</ul>

The two main R scripts are:

<ol>
          <li> `Modeling.R` - Fits the extremal dependence with provided data </li>
          <li> `Summary.R` - Provides a summary of the fitted model with related plots and tables  </li>
</ol>

Auxillary and utility functions are include in the R files:
<ol>
 <li>  `Simulate_data.R` - A simple simulation of max-stable processes with nonstationary extremal dependence;</li>
 <li>  `Algorithm1.R` - Function for Algorithm 1 described in the paper;</li>
 <li>  `Merge_subr.R` - Function for the subregion merging process described in the paper;</li>
 <li>  `Lambda_tuning.R` - Function for the lambda-tuning described in the paper;</li>
 <li> `Fit.R` - Some fitting function for convenience, using r-optim;</li>
 <li>  `Objectives.R` - Pairwise likelihood functions for the Brown-Resnick process (and inverted counterpart);</li>
 <li>  `Utils.R` - Various other utility functions.</li>
</ol>



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
