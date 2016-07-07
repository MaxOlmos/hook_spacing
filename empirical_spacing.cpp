
#include <TMB.hpp>
// lognormal density
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
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
  DATA_VECTOR(log_yobs); 	// the response variable
  DATA_INTEGER(Ngroup); // number of groups
  DATA_INTEGER(Nobs);	// number of observations
  DATA_VECTOR(day);	// day relative to start of fishing
  DATA_VECTOR(spacing);	// the independent variable spacing (ft)
  DATA_IVECTOR(group);	// index for group

  // Fixed effects
  PARAMETER(logcpue_mean);	   // Mean CPUE (asymptote)
  PARAMETER(logcpue_sd);	   // SD of mean CPUE
  PARAMETER(sigma_obs_mean);	   // mean Obs error, log scale
  PARAMETER(sigma_obs_sd);	   // sd of obs errors, log scale
  PARAMETER(gamma);		   // Impact of day, bounded (0,1)
  PARAMETER(beta);		   // non-linear effect of spacing
  // Random effects
  PARAMETER_VECTOR(logcpue);	  // site level random effects for logcpue
  PARAMETER_VECTOR(logsigma_obs); // observation error at each site

  // Objective function and bookkeeping
  using namespace density;
  // vector<Type> jnll_comp(3);	// joint -log-likelihood in components
  // jnll_comp.setZero();		// initialize at zero
  Type jnll=0;

  // Derived quantities
  vector<Type> sigma_obs(Ngroup);
  vector<Type> cpue(Ngroup);
  for(int i=0; i<Ngroup; i++){
    cpue(i)=exp(logcpue(i));
    sigma_obs(i)=exp(logsigma_obs(i));
  }

  // The model predictions for each observation
  vector<Type> ypred_i(Nobs);
  for( int i=0; i<Nobs; i++){
    ypred_i(i) = cpue(group(i))*exp(-day(i)*gamma)*(1-exp(-beta*spacing(i)));
  }

  // Probability of random effects on group
  for(int i=0; i<Ngroup; i++){
    jnll-= dnorm(logcpue_mean, logcpue(i), logcpue_sd, true);
    jnll-= dnorm(sigma_obs_mean, logsigma_obs(i), sigma_obs_sd, true);
  }
  // Probability of the data, given random effects (likelihood)
  for( int i=0; i<Nobs; i++){
    jnll-= dnorm(log(ypred_i(i)), log_yobs(i), sigma_obs(group(i)), true );
  }

  // Reporting
  vector<Type> spacing_pred(70);
  vector<Type> spacing_std(70);
  for(int i=0; i<70; i++){
    // predict spacing effect for day=0 and mean group level
    spacing_pred(i)= exp(logcpue_mean)*exp(-0)*(1-exp(-beta*(i+1)));
  }
  for(int i=0; i<70; i++){
    spacing_std(i)= spacing_pred(i)/spacing_pred(17);
  }
  ADREPORT(spacing_std);
  return jnll;
}
