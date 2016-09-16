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
  DATA_VECTOR(spacing_i);	// the independent variable spacing (ft)
  DATA_IVECTOR(site_i);		// index for site

  // Fixed effects
  PARAMETER(eta_mean);		// Mean density*q*alpha
  PARAMETER(eta_sd);		// SD of mean density
  PARAMETER(sigma_mean);	// mean Obs error
  PARAMETER(sigma_sd);		// sd of obs errors
  PARAMETER(gamma);		// Impact of day, bounded (0,1)
  PARAMETER(beta);		// non-linear effect of spacing
  PARAMETER(lambda);		// exponent that controls non-linearity = 1

  // Random effects
  PARAMETER_VECTOR(eta_s);	  // site level random effects for logcpue
  PARAMETER_VECTOR(sigma_s); // observation error at each site

  // Objective function and bookkeeping
  using namespace density;
  // vector<Type> jnll_comp(3);	// joint -log-likelihood in components
  // jnll_comp.setZero();		// initialize at zero
  Type jnll=0;


  // Predicted catches
  vector<Type> mu_i(n_i);
  for( int i=0; i<n_i; i++){
    mu_i(i) =
      // spacing
      (1-pow(exp(-beta*spacing_i(i)),lambda))*
      // hooks
      hooks_i(i)*
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
    Type test=10;
    jnll-= dlognorm(catch_i(i), log(mu_i(i)) , sigma_s(site_i(i)), true );
  }

  // Reporting
  // alpha tilde in paper
  Type max_ehook = 1/(1-pow(exp(-18*beta), lambda));
  vector<Type> hook_power(70);
  for(int i=0; i<70; i++){
    // predict spacing effect for day=0 and mean site level
    hook_power(i)=
      (1-pow(exp(-beta*(i+1)),lambda))/ (1-pow(exp(-beta*18),lambda));
  }
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
