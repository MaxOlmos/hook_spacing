//// Reanalysis of the Hamley and Skud (1978) data for hook spacing
//// experiments

#include <TMB.hpp>
// lognormal density
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs
  DATA_INTEGER(n_s);		// number of sites
  DATA_INTEGER(n_i);		// number of observations
  DATA_VECTOR(catch_i); 	// the response variable
  DATA_VECTOR(hooks_i);		// number of hooks
  DATA_VECTOR(day_i);		// day relative to start of fishing
  DATA_FACTOR(spacing_i);	// the independent variable spacing (ft)
  DATA_IVECTOR(site_i);		// index for site

  // Fixed effects
  PARAMETER(eta_mean);		// Mean density*q*alpha
  PARAMETER(eta_sd);		// SD of mean density
  PARAMETER(sigma_mean);	// mean Obs error
  PARAMETER(sigma_sd);		// sd of obs errors
  PARAMETER(gamma);			// Impact of day, bounded (0,1)
  PARAMETER(beta);			// non-linear effect of spacing
  PARAMETER(lambda);		// exponent that controls non-linearity = 1

  // Random effects
  PARAMETER_VECTOR(eta_s);	 // site level random effects for logcpue
  PARAMETER_VECTOR(sigma_s); // observation sd at each site

  // Objective function and bookkeeping
  using namespace density;
  Type jnll=0;

  // Calculate spacing effect from parameters
  int n_ft=45;
  vector<Type> spacing(n_ft);
  spacing.setZero();
  for(int ft=0; ft<n_ft; ft++){
    spacing(ft)=1-pow(exp(-beta*(ft+1)),lambda);
  }
  // Calculate relative hook power. Need to be careful with indexing
  // here. ft=0 => a 1 foot spacing. Thus hook_power(0) is hook power at 1
  // ft.
  vector<Type> hook_power(n_ft);
  for(int ft=0; ft<n_ft; ft++){
    hook_power(ft)=spacing(ft)/spacing(17);
  }

  // Predicted catches= effective hooks * density * q (q=1)
  vector<Type> mu_i(n_i);
  for(int i=0; i<n_i; i++){
    mu_i(i) =
      // Effective hooks
      hooks_i(i)* hook_power(spacing_i(i)-1)*
      // site level density
      exp(eta_s(site_i(i))-day_i(i)*gamma);
  }

  // Probability of random effects on site
  for(int s=0; s<n_s; s++){
    jnll-= dnorm(eta_mean, eta_s(s), eta_sd, true);
    jnll-= dnorm(sigma_mean, sigma_s(s), sigma_sd, true);
  }
  // Probability of the data, given random effects (likelihood)
  for( int i=0; i<n_i; i++){
    jnll-= dnorm(log(catch_i(i)), log(mu_i(i)) , sigma_s(site_i(i)), true );
  }

  // Reporting
  // max_ehook is f_infinity in paper
  Type max_ehook = 1/(1-pow(exp(-18*beta), lambda));
  REPORT(mu_i);
  ADREPORT(hook_power);
  ADREPORT(beta);
  ADREPORT(gamma);
  ADREPORT(lambda);
  ADREPORT(max_ehook);
  ADREPORT(sigma_mean);
  ADREPORT(sigma_sd);
  ADREPORT(eta_mean);
  ADREPORT(eta_sd);
  return jnll;
}
