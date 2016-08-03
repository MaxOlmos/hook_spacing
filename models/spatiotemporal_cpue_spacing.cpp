
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
  //// Data inputs
  DATA_INTEGER(likelihood); // the likelihood function to use
  DATA_INTEGER(form);	    // form of the hook spacing; 1=RE; 2=H&S
  DATA_INTEGER(space);	    // form of the spatial component; 0=NS, 1=S; 2=ST
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
  DATA_VECTOR(catch_i);  // catch/hook (response variable); natural scale
  DATA_VECTOR(hooks_i); // the number of hooks
  // SPDE objects from R-INLA
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  //// Parameters
  // Fixed effects
  PARAMETER(intercept);		   // global mean log catch
  PARAMETER_VECTOR(beta_year);	   // year effect
  PARAMETER_VECTOR(beta_geartype); // geartype effect
  PARAMETER_VECTOR(beta_month);	   // month effect
  PARAMETER_VECTOR(beta_hooksize); // hooksize effect
  PARAMETER(beta_depth);	   // depth effect
  PARAMETER(beta_spacing);	   // from H&S formula
  PARAMETER(lambda);		   // from H&S formula
  //PARAMETER(alpha_spacing);	   // from H&S formula
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
  int n_i = catch_i.size();	// number of observations
  Type nll_likelihood=0;	// likelihood of data
  Type nll_omega=0;		// spatial
  Type nll_epsilon=0;		// spatio-temporal
  Type nll_spacing=0;		// RW on hook spacing

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_O) * exp(2*ln_kappa));
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau_E) * exp(2*ln_kappa));
  Type Sigma = exp(ln_obs);
  Eigen::SparseMatrix<Type>
    Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  // Calculate the spacing effects by foot (random effects)
  vector<Type> spacing(n_ft);
  // ft here is offset by -1 so ft=0 => spacing of 1ft
  for(int ft=0; ft<n_ft; ft++){
    // Additive effect of hook spacing for a given foot
    if(form==1) {
      if(ft==0) spacing(ft)=spacing_devs(ft); // initialize at first dev
      else spacing(ft)=spacing(ft-1)+spacing_devs(ft);
    }
    // Multiplicative form used by H&S.
    if(form==2) spacing(ft)=(1-pow(exp(-beta_spacing*(ft+1)), lambda));
  }
  // Standardized effect of spacing
  vector<Type> spacing_std(n_ft);
  for(int ft=0; ft<n_ft; ft++)
    // ft is again offset by 1 here, so ft=0 =>  spacing of 1ft
    spacing_std(ft)=spacing(ft)/spacing(17);

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
  if(space>0) nll_omega += SCALE( GMRF(Q), 1/exp(ln_tau_O) )( omega_s );
  if(space>1) {
  for( int t=0; t<n_t; t++)	// spatio-temporal
    nll_epsilon += SCALE( GMRF(Q), 1/exp(ln_tau_E) )( epsilon_st.col(t) );
  }

  // hook spacing effects
  for(int ft=0; ft<n_ft; ft++)
    nll_spacing -= dnorm(Type(0.0), spacing_devs(ft),exp(ln_spacing), true);

  // Probability of the data, given random effects
  for( int i=0; i<n_i; i++){
    // Probability of data conditional on random effects (likelihood)
    if(likelihood==1) // lognormal case
      nll_likelihood -=
	dnorm(log(catch_i(i)), log(hooks_i(i))+mu_i(i), exp(ln_obs), true);
  }

  // Predict relative average catch rate over time
  vector<Type> cph_t(n_t);
  for( int t=0; t<n_t; t++)
    cph_t(t)=  exp(intercept + beta_year(t) + beta_month(0) +
		   beta_geartype(0) + beta_hooksize(0) + beta_depth*80);

  // Reporting
  Type jnll = nll_likelihood+nll_omega+nll_epsilon+nll_spacing;
  vector<Type> resids = log(hooks_i) +mu_i-log(catch_i);
  // REPORT(jnll_comp);
  REPORT(nll_likelihood);
  REPORT(nll_omega);
  REPORT(nll_epsilon);
  REPORT(nll_spacing);
  REPORT(intercept);
  REPORT(beta_depth);
  REPORT(spacing_std);
  REPORT(spacing);
  if(form==2){
    Type max_ehook = 1/(1-pow(exp(-18*beta_spacing), 1));
    ADREPORT(max_ehook);
    ADREPORT(beta_spacing);
    ADREPORT(lambda);
  }
  ADREPORT(Range);
  ADREPORT(SigmaE);
  ADREPORT(SigmaO);
  ADREPORT(Sigma);
  ADREPORT(cph_t);
  ADREPORT(spacing_std);
  REPORT(resids);
  return jnll;
}
