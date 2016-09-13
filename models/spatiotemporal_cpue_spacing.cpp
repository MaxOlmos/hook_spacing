
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
  DATA_INTEGER(n_s); // number of sites (grids)
  DATA_INTEGER(n_ft); // number of spacings.. 1:n_ft, with ft(0)=0 assumed
  DATA_INTEGER(n_v);  // number of vessels
  // Indices for factors
  DATA_FACTOR(spacing_i);
  DATA_FACTOR(year_i);
  DATA_FACTOR(geartype_i);
  DATA_FACTOR(month_i);
  DATA_FACTOR(hooksize_i);
  DATA_FACTOR(vessel_i);
  // DATA_FACTOR(statarea_i);
  // vectors of real data
  DATA_VECTOR(depth_i); // depth covariate
  DATA_VECTOR(catch_i);  // catch/hook (response variable); natural scale
  DATA_VECTOR(hooks_i); // the number of hooks
  DATA_VECTOR(area_s); // area for each grid cell, for area-weighted CPUE
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
  PARAMETER(beta_depth);	   // linear depth effect
  PARAMETER(beta_depth2);	   // quadratic depth effect
  PARAMETER(beta_spacing);	   // from H&S formula
  PARAMETER(alpha_spacing);	   // from H&S formula
  PARAMETER(lambda);		   // from H&S formula
  // Variances
  PARAMETER(ln_tau_O);		  // spatial process
  PARAMETER(ln_tau_E);		  // spatio-temporal process
  PARAMETER(ln_kappa);		  // decorrelation distance (kind of)
  PARAMETER(ln_obs);		  // Observation variance
  PARAMETER(ln_spacing);	  // Hook spacing effect
  PARAMETER(ln_vessel);		  // vessel effect
  // Random effects
  PARAMETER_VECTOR(spacing_devs); // hook spacing; n_ft length
  PARAMETER_VECTOR(vessel_v);	  // vessel effect; n_v length
  PARAMETER_VECTOR(omega_s);	  // spatial effects; n_s length
  // This is a matrix of (centers by years)
  PARAMETER_ARRAY(epsilon_st);	// spatio-temporal effects; n_s by n_t matrix

  // Objective function and bookkeeping
  using namespace density;
  int n_i = catch_i.size();	// number of observations
  Type nll_likelihood=0;	// likelihood of data
  Type nll_omega=0;		// spatial
  Type nll_epsilon=0;		// spatio-temporal
  Type nll_spacing=0;		// RW on hook spacing
  Type nll_vessel=0;		// vessel effect

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
    // smoother on effect of hook spacing for a given foot
    if(form==1) {
      if(ft==0) spacing(ft)=0; // initialize at first dev
      else spacing(ft)=spacing(ft-1)+spacing_devs(ft);
    }
    // Multiplicative form used by H&S.
    if(form==2) spacing(ft)=(1-pow(exp(-beta_spacing*(ft+1)), lambda));
    if(form==3) spacing(ft)=1; // no effect
  }
  // Calculate relative hook power
  vector<Type> hook_power(n_ft);
  for(int ft=0; ft<n_ft; ft++){
    // ft is again offset by 1 here, so ft=0 =>  spacing of 1ft
    if(form==1) hook_power(ft)=exp(spacing(ft))/exp(spacing(17));
    if(form==2) hook_power(ft)=spacing(ft)/spacing(17);
    if(form==3) hook_power(ft)=spacing(ft);
  }

  // The model predicted catch for each observation, in natural space:
  // catch=density*hook_power*hooks*catchability
  vector<Type> mu_i(n_i);
  for( int i=0; i<n_i; i++){
    mu_i(i) = hooks_i(i)*hook_power(spacing_i(i)-1)*
      exp(intercept + beta_year(year_i(i)) +
	  beta_month(month_i(i)) + beta_geartype(geartype_i(i)) +
	  beta_hooksize(hooksize_i(i)) +
	  beta_depth*depth_i(i) + beta_depth2*depth_i(i)*depth_i(i) +
	  vessel_v(vessel_i(i))+
	  omega_s(s_i(i)) + epsilon_st(s_i(i),year_i(i)));
  }

  //// Probability of random effects
  // Space and spatio-temporal
  if(space>0) nll_omega += SCALE(GMRF(Q), 1/exp(ln_tau_O))(omega_s);
  if(space>1) {
  for( int t=0; t<n_t; t++)
    nll_epsilon += SCALE(GMRF(Q), 1/exp(ln_tau_E))(epsilon_st.col(t));
  }
  // hook spacing effects
  for(int ft=0; ft<n_ft; ft++)
    nll_spacing -= dnorm(Type(0.0), spacing_devs(ft),exp(ln_spacing), true);
  // vessel effects
  for(int v=0; v<n_v; v++)
    nll_vessel -= dnorm(Type(0.0), vessel_v(v),exp(ln_vessel), true);

  //// Probability of the data, given random effects (likelihood)
  for( int i=0; i<n_i; i++){
    if(likelihood==1) // lognormal case; see top of file for function
      nll_likelihood -=
	dlognorm(catch_i(i), log(mu_i(i)), Sigma, true);
    if(likelihood==2) // gamma case
      nll_likelihood -=
    dgamma(catch_i(i), 1/pow(Sigma,2), mu_i(i)*pow(Sigma,2), true);
    if(likelihood==3) // inverse gaussian
      nll_likelihood-=
	dinvgauss(catch_i(i), mu_i(i), Sigma, true );
  }


  // Calculate joint negative log likelihood
  Type jnll = nll_likelihood+nll_omega+nll_epsilon+nll_spacing +nll_vessel;

  //// Derived quantities
  // catch per hook -- old way that doesn't make sense right now
  vector<Type> cph_t(n_t);
  for( int t=0; t<n_t; t++){
    cph_t(t)=  exp(intercept + beta_year(t) + beta_month(0) +
		   beta_geartype(0) + beta_hooksize(0) + beta_depth*80);
  }
  // Area weighted approach to estimating trends in abundance
  vector<Type> area_weighted_density_t(n_t); // one for each year
  for( int t=0; t<n_t; t++){
    area_weighted_density_t(t)=0;
    for(int s=0; s<n_s; s++){
      area_weighted_density_t(t)+=
	// add depth here??
	area_s(s)*exp(intercept+beta_year(t)+ omega_s(s) + epsilon_st(s,t));
    }
  }
  // Pearson residuals; (pred-obs)/sd(pred)
  vector<Type> resids;
  if(likelihood==1) resids=(log(mu_i)-log(catch_i))/Sigma;
  if(likelihood==2) resids=(mu_i-catch_i)/(Sigma*mu_i);
  if(likelihood==3) resids=(mu_i-catch_i)/(sqrt(mu_i+mu_i+mu_i)/sqrt(Sigma));
  vector<Type> preds = mu_i;

  //// Reporting
  REPORT(nll_likelihood);
  REPORT(nll_omega);
  REPORT(nll_epsilon);
  REPORT(nll_spacing);
  REPORT(nll_vessel);
  REPORT(intercept);
  REPORT(beta_depth);
  REPORT(hook_power);
  REPORT(spacing);
  REPORT(resids);
  REPORT(preds);
  if(form==2){
    Type max_ehook = 1/(1-pow(exp(-18*beta_spacing), lambda));
    ADREPORT(max_ehook);
    ADREPORT(beta_spacing);
    ADREPORT(alpha_spacing);
    ADREPORT(lambda);
  } else {
      ADREPORT(ln_spacing);	  // Hook spacing effect
  }
  ADREPORT(intercept);		   // global mean log catch
  ADREPORT(beta_year);	   // year effect
  ADREPORT(beta_geartype); // geartype effect
  // ADREPORT(beta_month);	   // month effect
  // ADREPORT(beta_hooksize); // hooksize effect
  ADREPORT(beta_depth);	   // linear depth effect
  ADREPORT(beta_depth2);   // quadratic depth effect
  ADREPORT(Sigma);	   	// observation
  // random effect variances
  Type Sigma_vessel=exp(ln_vessel);
  ADREPORT(Sigma_vessel);
  // Derived geospatial components
  ADREPORT(Range);		// geostatistical range
  ADREPORT(SigmaE);		// space
  ADREPORT(SigmaO);		// spatiotemporal
  // Spacing and CPUE calcs
  ADREPORT(hook_power);
  ADREPORT(cph_t);
  ADREPORT(area_weighted_density_t);
  return jnll;
}
