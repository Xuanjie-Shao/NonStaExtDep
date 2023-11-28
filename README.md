# Flexible Modeling of Nonstationary Extremal Dependence Using Spatially-Fused LASSO and Ridge Penalties
Functions for modeling nonstationary extremal dependence using locally-statioanry max-stable processes with LASSO and ridge penalization.
The provided code is in support of the paper:
<ul> 
          <li> Shao, X., Hazra, A., Richards, J., and Huser, R. (2023+). Flexible modeling of non-stationary extremal dependence using spatially-fused LASSO and ridge penalties. <u><a href="https://arxiv.org/abs/2210.05792" download>ArXiv</a></u> </li>
</ul>

The two main R scripts are:

<ol>
          <li> `Modeling.R` - Fits the extremal dependence with provided data </li>
          <li> `Summary.R` - Provides a summary of the fitted model with related plots and tables. See Figures 4-6 and Table 3 in the paper.  </li>
</ol>
which should be run sequentially. Running the code as is will perform model fitting and summary for a simulated data example. By changing the code `dat = "Simulated"` to `dat = "NepalExtended"` at the top of each script will run model fitting and summary for the data used in the application (see below). Note that the former takes only a few minutes to compile, whilst the latter will take a few hours!

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

Included in the repo are two Rdata files:
<ol>
 <li>  `NepalExtended.Rdata` - The gridded data of monthly maximum temperature dataset from Nepal and its surrounding Himalayan and sub-Himalayan regions used in the data application of the paper. This file includes the marginal parameter estiamtes (GEV parameters) derived using the Max-and-Smooth method: postman.mu (location), postman.sigma (scale), postman.xi (shape). </li>
 <li>  `Simulated.Rdata` - Simulated data from `Simulate_data.R`, including the coordinate and true (dependence) parameter information. The true partition is partition P1 mentioned in the paper. </li>
</ol>

Some further remarks:
<ul>
<li>  The current program only works for gridded data, but an extension to general lattice data is available.</li>
<li>  Nonstationarity is assumed for the input data.</li>
<li>  Input data must be renormalized to unit Fréchet margins to fit the max-stable processes.</li>
<li>  Some difficulties in range estimation may emerge.</li>
</ul> 

