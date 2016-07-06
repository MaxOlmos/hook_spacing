
#include <TMB.hpp>
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  // DATA_INTEGER( n_s );
  // DATA_INTEGER( n_t );

  //  DATA_VECTOR( a_s );  // Area associated with location s
  DATA_VECTOR( logcpue_i );  // counts for observation i
  DATA_FACTOR( s_i );  // Random effect index for observation i
  //  DATA_FACTOR( t_i );  // Random effect index for observation i
  DATA_FACTOR(year_i);
  DATA_FACTOR(geartype_i);
  DATA_FACTOR(month_i);
  DATA_FACTOR(hooksize_i);
  DATA_VECTOR(depth_i);

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER( intercept ); // global mean log cpue
  PARAMETER_VECTOR(beta_year); // year effect
  PARAMETER_VECTOR(beta_geartype); // geartype effect
  PARAMETER_VECTOR(beta_month); // month effect
  PARAMETER_VECTOR(beta_hooksize); // hooksize effect
  PARAMETER(beta_depth); // depth effect
  // Spatia variances
  PARAMETER( ln_tau_O );
  //  PARAMETER( ln_tau_E );
  PARAMETER( ln_kappa );
  // Observation variance
  PARAMETER( ln_obs);

  // Random effects
  PARAMETER_VECTOR( omega_s );
  // PARAMETER_ARRAY( epsilon_st );

  // Objective funcction
  using namespace density;
  int n_i = logcpue_i.size();
  //  int n_s = omega_s.size();
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_O) * exp(2*ln_kappa));
  //Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau_E) * exp(2*ln_kappa));

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  jnll_comp(1) += SCALE( GMRF(Q), 1/exp(ln_tau_O) )( omega_s );
  // for( int t=0; t<n_t; t++){
  //   jnll_comp(2) += SCALE( GMRF(Q), 1/exp(ln_tau_E) )( epsilon_st.col(t) );
  // }

  // True density and abundance
  //array<Type> log_d_st( n_s, n_t );
  vector<Type> mu_i(n_i);
  //  for( int t=0; t<n_t; t++){
  for( int i=0; i<n_i; i++){
    mu_i(i) = intercept +
      beta_year(year_i(i)) +
      beta_month(month_i(i)) +
      beta_geartype(geartype_i(i)) +
      beta_hooksize(hooksize_i(i)) +
      beta_depth*depth_i(i)+
      omega_s(s_i(i));// + epsilon_st(s,t);
    // Probability of data conditional on random effects
    if( !isNA(logcpue_i(i)) )
      jnll_comp(0) -= dnorm( logcpue_i(i), mu_i(i), exp(ln_obs), true );
  }

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( Range );
  //  REPORT( SigmaE );
  REPORT( SigmaO );
  REPORT( mu_i );

  return jnll;
}
