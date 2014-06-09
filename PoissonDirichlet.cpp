/*    PoissonDirichlet:
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *                       Timothy Daley
 *
 *    Authors: Andrew D. Smith and Timothy Daley
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "PoissonDirichlet.hpp"

#include <smithlab_utils.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_rng.h>

#include <fstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numeric>


using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;

using std::max;
using std::fabs;
using std::accumulate;



//////////////////////////////////////////////////////////////////////
// 2 parameter Poisson Dirichlet Distribution

const double PDD::tolerance = 1e-20;
const size_t PDD::max_iter = 50;


static double
empirical_bayes_log_eq(const vector<double> &counts,
		       const double counts_sum,
		       const double theta,
		       const double sigma){
  double log_return_val = 0.0;

  for(size_t i = 1; i < counts.size(); i++)
    log_return_val += log(theta + i*sigma);

  log_return_val -= gsl_sf_lngamma(theta + 1.0);

  log_return_val += gsl_sf_lngamma(theta + counts_sum);

  for(size_t i = 0; i < counts.size(); i++)
    log_return_val += gsl_sf_lngamma(counts[i] - sigma) - gsl_sf_lngamma(1.0 - sigma);

  return log_return_val;
}


static double
empirical_bayes_log_dsigma(const vector<double> &counts,
			   const double counts_sum,
			   const double theta,
			   const double sigma){
  double log_return_val = 0.0;

  for(size_t i = 1; i < counts.size(); i++)
    log_return_val += static_cast<double>(i)/(theta + i*sigma);

  for(size_t i = 0; i < counts.size(); i++)
    log_return_val -= (gsl_sf_psi(counts[i] - sigma) - gsl_sf_psi(1.0 - sigma));

  return log_return_val;
}

static double
empirical_bayes_log_dtheta(const vector<double> &counts,
			   const double counts_sum,
			   const double theta,
			   const double sigma){
  double log_return_val = 0.0;

  for(size_t i = 1; i < counts.size(); i++)
    log_return_val += 1.0/(theta + i*sigma);

  log_return_val += gsl_sf_psi(theta + 1.0);

  log_return_val -= gsl_sf_psi(theta + counts_sum);

  return log_return_val;
}

static double
empirical_bayes_log_dtheta_dsigma(const vector<double> &counts,
				  const double counts_sum,
				  const double theta,
				  const double sigma){
  double log_return_val = 0.0;

  for(size_t i = 1; i < counts.size(); i++)
    log_return_val -= static_cast<double>(i)/pow(theta + i*sigma, 2);

  return log_return_val;
}

static double
empirical_bayes_log_d2theta(const vector<double> &counts,
				  const double counts_sum,
				  const double theta,
				  const double sigma){
  double log_return_val = 0.0;

  for(size_t i = 1; i < counts.size(); i++)
    log_return_val -= 1.0/pow(theta + i*sigma, 2);

  log_return_val += gsl_sf_psi_1(theta + 1.0);

  log_return_val -= gsl_sf_psi_1(theta + counts_sum);

  return log_return_val;
}

static double
empirical_bayes_log_d2sigma(const vector<double> &counts,
				  const double counts_sum,
				  const double theta,
				  const double sigma){
  double log_return_val = 0.0;

  for(size_t i = 1; i < counts.size(); i++)
    log_return_val -= pow(i/(theta + i*sigma), 2);

  for(size_t i = 0; i < counts.size(); i++)
    log_return_val += gsl_sf_psi_1(counts[i] - sigma) - gsl_sf_psi_1(1.0 - sigma);

  return log_return_val;
}

static double
square_matrix_determinant(const vector< vector<double> > &mat){
  return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
}

static void
invert_square_matrix(vector< vector<double> > &mat){
  assert(mat.size() == 2 && mat[0].size() == 2);
  const vector< vector<double> > orig_mat(mat);
  const double det = square_matrix_determinant(orig_mat);
  mat[0][0] = orig_mat[1][1]/det;
  mat[0][1] = -orig_mat[0][1]/det;
  mat[1][0] = -orig_mat[1][0]/det;
  mat[1][1] = orig_mat[0][0]/det;
}

static void
NewtonRaphson2dStep(const vector<double> &counts,
		    const double counts_sum,
		    const double current_theta,
		    const double current_sigma,
		    vector<double> &step){
  // set derivative
  vector<double> deriv(2, 0.0);
  deriv[0] = empirical_bayes_log_dtheta(counts, counts_sum,
					current_theta, current_sigma);
  deriv[1] = empirical_bayes_log_dsigma(counts, counts_sum,
					current_theta, current_sigma);
  // set jacobian
  vector< vector<double> > jacobian(2, vector<double>(2, 0.0));
  jacobian[0][0] = 
    empirical_bayes_log_d2theta(counts, counts_sum, current_theta, current_sigma);
  jacobian[0][1] = 
    empirical_bayes_log_dtheta_dsigma(counts, counts_sum, current_theta, current_sigma);
  jacobian[1][0] = 
    empirical_bayes_log_dtheta_dsigma(counts, counts_sum, current_theta, current_sigma);
  jacobian[1][1] = 
    empirical_bayes_log_d2sigma(counts, counts_sum, current_theta, current_sigma);

  // jacobian is now the inverse jacobian
  invert_square_matrix(jacobian);

  // update step
  step.clear();
  step.push_back(deriv[0]*jacobian[0][0] + deriv[1]*jacobian[1][0]);
  step.push_back(deriv[0]*jacobian[0][1] + deriv[1]*jacobian[1][1]);
}

static inline bool
check_sigma(const double sigma){
  if(sigma >= 0.0 && sigma <= 1.0)
    return true;
  //else
  return false;
}

static inline bool
check_theta(const double theta){
  if(theta >= 0.0)
    return true;
  //else
  return false;
}

void
PDD::newton_raphson_estim_params(const bool VERBOSE, const vector<double> &counts){
  const double counts_sum = 
    accumulate(counts.begin(), counts.end(), 0.0);

  // start at random point
   // Setup the random number generator
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default); 
  srand(time(0) + getpid());
  gsl_rng_set(rng, rand());
  double prev_sigma = sigma;
  double prev_theta = theta;

  double current_theta = theta;
  double current_sigma = sigma;

  double error = std::numeric_limits<double>::max();
  size_t iter = 0;
  do {
    prev_theta = current_theta;
    prev_sigma = current_sigma;
    vector<double> movement;
    NewtonRaphson2dStep(counts, counts_sum, prev_theta, prev_sigma, movement);
    current_theta = prev_theta - movement[0];
    current_sigma = prev_sigma - movement[1];

    // if we go out of bounds, modify the step size
    if(!(check_theta(current_theta) && check_sigma(current_sigma))){
      double lower_multiplier = std::min(fabs(prev_theta/movement[0]), fabs(prev_sigma/movement[1]))/2.0;
      double upper_multiplier = (current_sigma > 1.0) ? ((prev_sigma - 1.0)/movement[1]) : std::numeric_limits<double>::max();
      double multiplier = std::min(lower_multiplier, upper_multiplier);

      current_theta = prev_theta - movement[0]*multiplier;
      current_sigma = prev_sigma - movement[1]*multiplier;

      if(VERBOSE){
	cerr << "multiplier = " << multiplier << endl;
	cerr << "step     = (" << movement[0] << ", " << movement[1] << ")" << endl;
	cerr << "new step = (" << movement[0]*multiplier << ", " << movement[1]*multiplier << ")" << endl; 
	cerr << "params     = (" << prev_theta << ", " << prev_sigma << ")" << endl;
	cerr << "new params = (" << current_theta << ", " << current_sigma << ")" << endl;
      }
    }

    error = 
      pow(empirical_bayes_log_eq(counts, counts_sum, current_theta, current_sigma)
	  - empirical_bayes_log_eq(counts, counts_sum, prev_theta, prev_sigma), 2);
    iter++;

    if(VERBOSE){
      cerr << "iter = " << iter << endl;
      cerr << "func = " << empirical_bayes_log_eq(counts, counts_sum, current_theta, current_sigma) << endl;
      cerr << "(theta, sigma) = (" << current_theta << ", " << current_sigma << ")" << endl;
    }
  } while (iter < max_iter && error > tolerance);

  sigma = current_sigma;
  theta = current_theta;

  if(sigma < min_allowed_sigma)
    cerr <<  "sigma = " << sigma << endl;

  if(sigma > max_allowed_sigma)
    cerr << "sigma = " << sigma << endl;

  if(theta < min_allowed_theta)
    cerr << "theta = " << theta << endl;

  assert(sigma >= min_allowed_sigma && sigma <= max_allowed_sigma && theta >= min_allowed_theta);
}

// see equation 6 of Favaro et al. J.R.Stat.Soc (2009)
double
PDD::expected_additional_distinct(const size_t current_total_count,
				  const size_t current_distinct,
				  const size_t future_total_count){
  const double ascending_factorial_term = 
    exp(gsl_sf_lngamma(theta + current_total_count + sigma + future_total_count)
	- gsl_sf_lngamma(theta + current_total_count + sigma)
	- gsl_sf_lngamma(theta + current_total_count + future_total_count)
	+ gsl_sf_lngamma(theta + current_total_count));

  return (current_distinct + theta/sigma)*(ascending_factorial_term - 1.0);
}

// see equation 7 of Favaro et al. J.R.Stat.Soc (2009)
double 
PDD::expected_discovery_prob(const size_t current_total_count,
			     const size_t current_distinct,
			     const size_t future_total_count){
  const double log_ascending_factorial_term = 
    gsl_sf_lngamma(theta + current_total_count + sigma + future_total_count)
    - gsl_sf_lngamma(theta + current_total_count + sigma)
    - gsl_sf_lngamma(theta + current_total_count + 1 + future_total_count)
    + gsl_sf_lngamma(theta + current_total_count + 1);

  return exp(log(theta + current_distinct*sigma) 
	     - log(theta + current_total_count) 
	     + log_ascending_factorial_term);
}

void
PDD::sample_PoissDir_counts(const gsl_rng *rng,
			    const size_t total_count,
			    vector<double> &sampled_counts){
  assert(total_count > 0);
  sampled_counts.clear();

  sampled_counts.push_back(1.0);
  size_t current_total_counts = 1;
  while(current_total_counts < total_count){
    double u = gsl_rng_uniform(rng);
    double test_val = 0.0;
    // add a count to x[i] with probability (x[i] - sigma)/(theta + current_total_counts)
    for(size_t i = 0; i < sampled_counts.size(); i++){
      test_val += (sampled_counts[i] - sigma)/(theta + current_total_counts);
      if (u < test_val){
	sampled_counts[i]++;
	break;
      }
    }
    // add new count with probability (theta + sampled_counts.size())/(theta + current_total_counts)
    if (u > test_val){
      sampled_counts.push_back(1.0);
    }

    current_total_counts++;
  }
}
