These scripts provide non-linear exponential decay fits of immune protection.

**Nonlinear_Meta_Regression_Script.R** runs the full meta-regression analysis and calculates the Bayesian Information Criterion (BIC) associated with each method.
Estimated parameters, BIC, log-likelihoods, etc are stored as "Non_Linear_fit.rds". 

**Filter_Parameters.R** Takes as input the a "Non_Linear_fit.rds" file and selects the best fitting models incorporating BIC plus some out of sample fitting.
Fitting statistics are stored as "Fitting_Summary.rds" and a reduced parameter table as "Non_Linear_Fit_reduced.rds"

**Simplified_Nonlinear_Meta_Regression_Script.R** runs a single model for each outcome. The models are customized to give parsimonious and biologically plausible outputs. These are the parameters which we intend to use in our modeling. The outputs are stored as "Non_Linear_Fit_simple.rds".