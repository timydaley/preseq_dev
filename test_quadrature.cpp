/*    test_quadrature:
 *
 *    Copyright (C) 2014 University of Southern California and
 *			 Andrew D. Smith and Timothy Daley
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
 *    along with this program.	If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <smithlab_os.hpp>

#include "moment_sequence.hpp"
#include "ZTNB.hpp"
#include "library_size_estimates.hpp"
#include "newtons_method.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::fixed;
using std::setprecision;

/*
void
generate_NBD(const vector<double> &mus, 
	     const vector<double> &alphas,
	     const vector<double> &mixing,
	     vector<size_t> &sample){

  const gsl_rng_type *T;
  gsl_rng *rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL) + getpid());

  for(size_t i = 0; i < sample.size(); i++){
    double u = gsl_rng_uniform(rng);
    double test_val = 0.0;
    for(size_t j = 0; j < mixing.size(); j++){
      test_val += mixing[j];
      if(u <= test_val){
	const double n = 1/alphas[j];
	const double p = n/(n + mus[j]);
	sample[i] = gsl_ran_negative_binomial(rng, p, n);
	break;
      }
    }
  }
     
} 
*/


/*
// m[j] = 'th modified moment, v[j]=j'th moment
// monic generalized laguerre polynomial w/ k = 1/alpha,
// p = mu*alpha/(1+ mu*alpha) : l_{j}(x)
// orthogonal to x e^{-x}
// l_j(x) = \sum_l=0^j j!/l! binom{1+j}{j-l} (-1)^{l+j} x^l
// m[j] = \sum_{l=0}^j j!/l! binom{1+j}{j-l} (-1)^{l+j} v[l]
static void
laguerre_modified_moments(const vector<double> &orig_moments,
			  const double mu,
			  const double alpha,
			  const size_t n_points,
			  vector<double> &modified_moments){
  modified_moments.resize(2*n_points, 0.0);
  const double k = 1/alpha;
  const double p = mu*alpha/(1.0 + mu*alpha);
  for(size_t n = 0; n < modified_moments.size(); n++){
    for(size_t l = 0; l <= n; l++){
      modified_moments[n] += 
	exp(gsl_sf_lngamma(n + k + 1) - gsl_sf_lngamma(n - l + 1)
	    - gsl_sf_lngamma(k + l) + gsl_sf_lnfact(n)
	    - gsl_sf_lnfact(l + 1) + (n - l)*log(p)
	    + l *log(orig_moments[l]))*pow(-1, n + l); 
    } 
  } 
}
*/


void
expected_NegBin_counts_hist(const double mu,
			    const double alpha,
			    const double tolerance,
			    const size_t lib_size,
			    vector<double> &expected_counts_hist){
  expected_counts_hist.clear();
  expected_counts_hist.push_back(0.0);
  double count = 1;
  do{
    double expected_count = 
      exp(log(lib_size) + gsl_sf_lngamma(count + alpha) 
	  - gsl_sf_lngamma(count + 1.0) - gsl_sf_lngamma(alpha) 
	  + count*(log(alpha) + log(mu) - log(1.0 + alpha*mu))
	  - (log(1.0 + alpha*mu))/alpha);
    expected_counts_hist.push_back(expected_count);
    ++count;
  } while (expected_counts_hist.back() > tolerance);
}

int
main(const int argc, const char **argv) {

  try {

    /* FILES */
    string outfile;
    
    size_t num_points = 100;
    size_t lib_size = 1000000;
    double tolerance = 1e-20;
    size_t max_iter = 1000;
    size_t hist_max_terms = 1000;
    size_t bootstraps = 1000;
    double CI = 0.95;
    double mu = 1.0;
    double alpha = 1.0;
    double lower_prob = 1e-10;

    
    /* FLAGS */
    bool VERBOSE = false;
    //	bool SMOOTH_HISTOGRAM = false;	
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "",
			   "<sorted-bed-file>");
    opt_parse.add_opt("outfile", 'o', "output file for estimates",
		      false , outfile);
    //  opt_parse.add_opt("n_points",'p', "number of points for approximation",
    //		      false, num_points);
    opt_parse.add_opt("hist_max_terms",'h',"max terms in histogram",
		      false, hist_max_terms);
    opt_parse.add_opt("lib_size",'L', "library size",
		      false, lib_size);
    opt_parse.add_opt("mean", 'm', "mean of Negative Binomial",
		      false, mu);
    opt_parse.add_opt("alpha",'a',"dispersion parameter for Negative Binomial",
		      false, alpha);
    opt_parse.add_opt("tol",'t',"numerical tolerance",
		      false, tolerance);
    opt_parse.add_opt("max_iter",'i',"maximum # iterations",
		      false, max_iter);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps to perform",
		      false, bootstraps);
    opt_parse.add_opt("lower_prob", 'l', "lower bound on probability",
		      false, lower_prob);
    //	opt_parse.add_opt("terms",'t',"maximum number of terms", false, 
    //	     orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false, VERBOSE);
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    /******************************************************************/

    vector<double> expected_counts_hist;
    expected_NegBin_counts_hist(mu, alpha, 0.01, lib_size, expected_counts_hist);

    if(VERBOSE){
      cerr << "EXPECTED COUNTS = " << endl;
      for(size_t i = 0; i < expected_counts_hist.size(); i++)
	cerr << i << '\t' << expected_counts_hist[i]  << endl;
      cerr << endl;

      cerr << "EXPECTED_UNOBSERVED = "
	   <<  exp(log(lib_size) - (log(1.0 + alpha*mu))/alpha)
	   << endl;
    }

    cerr << "CHAO_LOWER_BOUND = " 
	 << chao87_lowerbound_unobserved(expected_counts_hist) << endl;

    cerr << "HARRIS_3MOMENT_LOWER_BOUND = "
	 << harris_3moments_unobserved(VERBOSE, expected_counts_hist) << endl;

    size_t n_points = 2;
    cerr << "2_POINT_QUADRATURE_LOWER_BOUND = "
	 << quadrature_unobserved_lower_bound(VERBOSE, expected_counts_hist,
					      tolerance, max_iter, n_points)
	 << endl;

    double lower_lambda = mu*lower_prob;
    cerr << "lower_limit = " << lower_lambda << endl;
    /*
    cerr << "2_POINT_QUADRATURE_UPPER_BOUND = "
	 << quadrature_unobserved_upper_bound(VERBOSE, expected_counts_hist,
					      tolerance, max_iter, 
					      lower_lambda, n_points)
	 << endl;
    */
    n_points = 3;
    cerr << "3_POINT_QUADRATURE_LOWER_BOUND = "
	 << quadrature_unobserved_lower_bound(VERBOSE, expected_counts_hist,
					      tolerance, max_iter, n_points)
	 << endl;
  
    n_points = 3;
    cerr << "3_POINT_QUADRATURE_UPPER_BOUND = "
	 << quadrature_unobserved_upper_bound(VERBOSE, expected_counts_hist,
					      tolerance, max_iter, 
					      lower_lambda, n_points)
	 << endl;

    n_points = 4;
    cerr << "4_POINT_QUADRATURE_LOWER_BOUND = "
	 << quadrature_unobserved_lower_bound(VERBOSE, expected_counts_hist,
					      tolerance, max_iter, n_points)
	 << endl;

    n_points = 4;
    cerr << "4_POINT_QUADRATURE_UPPER_BOUND = "
	 << quadrature_unobserved_upper_bound(VERBOSE, expected_counts_hist,
					      tolerance, max_iter, 
					      lower_lambda, n_points)
	 << endl;

    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////

  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
