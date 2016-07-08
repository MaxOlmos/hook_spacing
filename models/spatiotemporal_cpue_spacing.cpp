
#include <TMB.hpp>
// lognormal density
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}
// inverse gamma
template<class Type>
Type dinvgauss(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.5*log(shape) - 0.5*log(2*M_PI*pow(x,3)) - (shape * pow(x-mean,2) / (2*pow(mean,2)*x));
  if(give_log) return logres; else return exp(logres);
}
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Model code
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs
  DATA_INTEGER(likelihood); // the likelihood function to use
  DATA_INTEGER(form);	    // form of the hook spacing; 1=RE; 2=H&S
  DATA_FACTOR(s_i);  // Random effect index for observation i
  DATA_INTEGER(n_t); // number of years
  DATA_INTEGER(n_ft); // number of spacings.. 1:n_ft, with ft(0)=0 assumed
  // Indices for factors
  DATA_FACTOR(spacing_i);
  DATA_FACTOR(year_i);
  DATA_FACTOR(geartype_i);
  DATA_FACTOR(month_i);
  DATA_FACTOR(hooksize_i);
  // DATA_FACTOR(statarea_i);
  // vectors of real data
  DATA_VECTOR(depth_i); // depth covariate
  DATA_VECTOR(cph_i);  // catch/hook (response variable); natural scale
  // SPDE objects from R-INLA
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  // Fixed effects
  PARAMETER(intercept);		   // global mean log cph
  PARAMETER_VECTOR(beta_year);	   // year effect
  PARAMETER_VECTOR(beta_geartype); // geartype effect
  PARAMETER_VECTOR(beta_month);	   // month effect
  PARAMETER_VECTOR(beta_hooksize); // hooksize effect
  PARAMETER(beta_depth);	   // depth effect
  PARAMETER(beta_spacing);	   // from H&S formula
  PARAMETER(alpha_spacing);	   // from H&S formula
  // Variances
  PARAMETER(ln_tau_O);		  // spatial process
  PARAMETER(ln_tau_E);		  // spatio-temporal process
  PARAMETER(ln_kappa);		  // decorrelation distance (kind of)
  PARAMETER(ln_obs);		  // Observation variance
  PARAMETER(ln_spacing);	  // Hook spacing effect
  // Random effects
  PARAMETER_VECTOR(spacing_devs); // hook spacing; n_ft length
  PARAMETER_VECTOR(omega_s);	  // spatial effects
  // This is a matrix of (centers by years)
  PARAMETER_ARRAY(epsilon_st);	// spatio-temporal effects

  // Objective function and bookkeeping
  using namespace density;
  int n_i = cph_i.size();	// number of observations
  vector<Type> jnll_comp(3);	// joint -log-likelihood in components
  jnll_comp.setZero();		// initialize at zero

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_O) * exp(2*ln_kappa));
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau_E) * exp(2*ln_kappa));
  Type Sigma = exp(ln_obs);
  Eigen::SparseMatrix<Type>
    Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  // Calculate the spacing effects by foot (random effects)
  vector<Type> spacing(n_ft+1); // n+1 since ft=0 is included
  spacing(0)=0; // this is for ft=0
  for(int ft=1; ft<=n_ft; ft++){
    // Additive effect of hook spacing for a given foot
    if(form==1) spacing(ft)=spacing(ft-1)+spacing_devs(ft-1);
    // Multiplicative form used by H&S
    if(form==2) spacing(ft)=alpha_spacing*(1-exp(-beta_spacing*(ft)));
  }
  // Standardized effect of spacing
  vector<Type> spacing_std(n_ft+1);
  for(int ft=0; ft<=n_ft; ft++)
    spacing_std(ft)=spacing(ft)/spacing(18);

  // The model predictions for each observation
  vector<Type> mu_i(n_i);
  for( int i=0; i<n_i; i++){
    if(form==1)
    mu_i(i) = intercept +
      spacing(spacing_i(i))+ beta_year(year_i(i)) +
      beta_month(month_i(i)) + beta_geartype(geartype_i(i)) +
      beta_hooksize(hooksize_i(i)) + beta_depth*depth_i(i) +
      omega_s(s_i(i)) + epsilon_st(s_i(i),year_i(i));
    if(form==2)
      mu_i(i) = spacing(spacing_i(i))*
	(intercept + beta_year(year_i(i)) +
	 beta_month(month_i(i)) + beta_geartype(geartype_i(i)) +
	 beta_hooksize(hooksize_i(i)) + beta_depth*depth_i(i) +
	 omega_s(s_i(i)) + epsilon_st(s_i(i),year_i(i)));
  }

  // Probability of random effects
  jnll_comp(1) += // space
    SCALE( GMRF(Q), 1/exp(ln_tau_O) )( omega_s );
  for( int t=0; t<n_t; t++)	// spatio-temporal
    jnll_comp(2) += SCALE( GMRF(Q), 1/exp(ln_tau_E) )( epsilon_st.col(t) );
    // hook spacing effects
  Type test=0;		// for some reason had to do as separate variable??
  for(int ft=0; ft<n_ft; ft++)
    test -= dnorm(Type(0.0), spacing_devs(ft),exp(ln_spacing), true);

  // Probability of the data, given random effects
  for( int i=0; i<n_i; i++){
    // Probability of data conditional on random effects (likelihood)
    if(likelihood==1) // lognormal case
      jnll_comp(0) -= dlognorm(cph_i(i), mu_i(i), exp(ln_obs), true );
    else if(likelihood==2) // gamma case
      jnll_comp(0) -= dgamma(cph_i(i), mu_i(i), exp(ln_obs), true );
    else error("bad likelihood input");
  }

  // Predict relative average catch rate over time
  vector<Type> cph_t(n_t);
  for( int t=0; t<n_t; t++)
    cph_t(t)=  exp(intercept + beta_year(t) + beta_month(0) +
		   beta_geartype(0) + beta_hooksize(0) + beta_depth*80);

  // Reporting
  Type jnll = jnll_comp.sum() + test;
  vector<Type> resids = mu_i-cph_i;
  // REPORT(jnll_comp);
  REPORT(jnll);
  REPORT(Range);
  REPORT(SigmaE);
  REPORT(SigmaO);
  REPORT(Sigma);
  REPORT(intercept);
  REPORT(beta_depth);
  REPORT(spacing_std);
  REPORT(spacing);
  ADREPORT(cph_t);
  ADREPORT(spacing_std);
  REPORT(resids);
  return jnll;
}
