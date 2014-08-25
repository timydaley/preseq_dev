
/*    test_harris:
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
#include <gsl/gsl_cdf.h>


#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <smithlab_os.hpp>

#include "moment_sequence.hpp"
#include "ZTNB.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::fixed;
using std::setprecision;
using std::isfinite;
using std::log;


void
generate_NBD(const double mu,
	     const double alpha,
	     vector<size_t> &sample){

  const gsl_rng_type *T;
  gsl_rng *rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL) + getpid());
  const double n = 1/alpha;
  const double p = 1.0/(1.0 + mu*alpha);

  for(size_t i = 0; i < sample.size(); i++)
    sample[i] = gsl_ran_negative_binomial(rng, p, n);
     
} 


/*
// check 3 term recurrence to avoid non-positive elements
// truncate if non-positive element found
static void
check_three_term_relation(vector<double> &a,
			  vector<double> &b){

  // first entry is zero! Abort
  if(a[0] <= 0.0){
    a.clear();
    b.clear();
  }

  for(size_t i = 0; i < b.size(); i++){
    if(b[i] <= 0.0 || !isfinite(b[i])
       || a[i + 1] <= 0.0 || !isfinite(a[i + 1])){
      b.resize(i);
      a.resize(i + 1);
      break;
    }
  }
}
*/


// vals_hist[j] = n_{j} = # (counts = j)
// vals_hist_distinct_counts[k] = kth index j s.t. vals_hist[j] > 0
// stores kth index of vals_hist that is positive
// distinct_counts_hist[k] = vals_hist[vals_hist_distinct_counts[k]]
// stores the kth positive value of vals_hist
void
resample_hist(const gsl_rng *rng, const vector<size_t> &vals_hist_distinct_counts,
              const vector<double> &distinct_counts_hist,
              vector<double> &out_hist) {

  vector<unsigned int> sample_distinct_counts_hist(distinct_counts_hist.size(), 0);

  const unsigned int distinct =
    static_cast<unsigned int>(accumulate(distinct_counts_hist.begin(),
                                         distinct_counts_hist.end(), 0.0));

  gsl_ran_multinomial(rng, distinct_counts_hist.size(), distinct,
                      &distinct_counts_hist.front(),
                      &sample_distinct_counts_hist.front());

  out_hist.clear();
  out_hist.resize(vals_hist_distinct_counts.back() + 1, 0.0);
  for(size_t i = 0; i < sample_distinct_counts_hist.size(); i++)
    out_hist[vals_hist_distinct_counts[i]] =
      static_cast<double>(sample_distinct_counts_hist[i]);
}


void
quadrature_bootstraps(const bool VERBOSE,
		      const vector<double> &orig_hist,
		      const size_t bootstraps,
		      const size_t max_num_points,
		      const double tolerance,
		      vector<double> &quad_estimates){
  quad_estimates.clear();

  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand());


  vector<size_t> orig_hist_distinct_counts;
  vector<double> distinct_orig_hist;
  for (size_t i = 0; i < orig_hist.size(); i++){
    if (orig_hist[i] > 0) {
      orig_hist_distinct_counts.push_back(i);
      distinct_orig_hist.push_back(orig_hist[i]);
    }
  }

  const size_t max_iter = 10*bootstraps;
  for(size_t iter = 0; iter < max_iter && quad_estimates.size() < bootstraps; ++iter){

    vector<double> hist;
    resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);

  // initialize moments, 0th moment is 1
    vector<double> bootstrap_moments(1, 1.0);
  // mu_r = (r + 1)! n_{r+1} / n_1
    for(size_t i = 0; i < 2*max_num_points; i++)
      bootstrap_moments.push_back(exp(gsl_sf_lnfact(i + 2) 
				      + log(hist[i + 2])
				      - log(hist[1])) );

    MomentSequence bootstrap_mom_seq(bootstrap_moments);

    vector<double> bootstrap_points, bootstrap_weights;
    bool QUAD_SUCCESS = 
      bootstrap_mom_seq.QR_quadrature_rules(VERBOSE, max_num_points, tolerance,
					    max_iter, bootstrap_points,
					    bootstrap_weights);

    if(QUAD_SUCCESS && bootstrap_points.size() == max_num_points){

      const double weights_sum = accumulate(bootstrap_weights.begin(), 
					    bootstrap_weights.end(), 0.0);
      if(weights_sum != 1.0){
	for(size_t i = 0; i < bootstrap_weights.size(); i++)
	  bootstrap_weights[i] = bootstrap_weights[i]/weights_sum;
      }

      double estimated_integral = 0.0;
      for(size_t i = 0; i < bootstrap_weights.size(); i++)
	estimated_integral += hist[1]*bootstrap_weights[i]/bootstrap_points[i];

      quad_estimates.push_back(estimated_integral);
    }

  }

}

void
log_mean_quad(const bool VERBOSE,
	      const vector<double> &quad_estimates,
	      const double ci_level,
	      double &log_mean, 
	      double &log_lower_ci,
	      double &log_upper_ci){
  vector<double> log_quad_estimates(quad_estimates);
  for(size_t i = 0; i < log_quad_estimates.size(); i++)
    log_quad_estimates[i] = log(log_quad_estimates[i]);

  log_mean = exp(gsl_stats_mean(&log_quad_estimates[0], 1,
				log_quad_estimates.size()) );

  double log_std_dev = std::sqrt(gsl_stats_variance(&log_quad_estimates[0], 1, 
						    log_quad_estimates.size()) );

  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv((1.0 - ci_level)/2.0);
  log_lower_ci = exp(log(log_mean) - inv_norm_alpha*log_std_dev);
  log_upper_ci = exp(log(log_mean) + inv_norm_alpha*log_std_dev);
}



int
main(const int argc, const char **argv) {

  try {

    /* FILES */
    string quad_outfile;
    
    size_t num_points = 100;
    size_t lib_size = 1000000;
    double tolerance = 1e-20;
    size_t max_iter = 1000;
    size_t hist_max_terms = 1000;
    size_t bootstraps = 100;
    double ci_level = 0.95;
    double distro_alpha = 1.0;
    double distro_mu = 1.0;

    
    /* FLAGS */
    bool VERBOSE = false;
    //	bool SMOOTH_HISTOGRAM = false;	
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "",
			   "<sorted-bed-file>");
    opt_parse.add_opt("quad_outfile", 'q', "output file for quadrature estimates",
		      false, quad_outfile);
    opt_parse.add_opt("n_points",'p', "number of points for approximation",
		      false, num_points);
    opt_parse.add_opt("hist_max_terms",'h',"max terms in histogram",
		      false, hist_max_terms);
    opt_parse.add_opt("lib_size",'l', "library size",
		      false, lib_size);
    opt_parse.add_opt("mean", 'm', "mu for NegBin dist", false, distro_mu);
    opt_parse.add_opt("alpha",'a',"alpha for NegBin dist",
    		      false, distro_alpha);
    opt_parse.add_opt("tol",'t',"numerical tolerance",
    		      false, tolerance);
    opt_parse.add_opt("max_iter",'i',"maximum # iterations",
		      false, max_iter);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps to perform",
		      false, bootstraps);
    opt_parse.add_opt("ci_level",'c', "Confidence level",
		      false, ci_level);
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


    // BUILD THE HISTOGRAM
    //    double mu = sampled_reads/lib_size;
      if(VERBOSE)
	cerr << "GENERATE SAMPLE" << endl; 
      vector<size_t> sample_counts(lib_size, 0);
      generate_NBD(distro_mu, distro_alpha, sample_counts);
      const size_t max_observed_count = *std::max_element(sample_counts.begin(), sample_counts.end());
      vector<double> counts_hist(max_observed_count + 1, 0.0);
      for(size_t i = 0; i < sample_counts.size(); i++)
	counts_hist[sample_counts[i]]++;

      counts_hist[0] = 0;

      const double distinct_reads = accumulate(counts_hist.begin(), counts_hist.end(), 0.0);

      if (VERBOSE) {
	cerr << "LIBRARY_SIZE = " << lib_size << endl;
	cerr << "MU = " << distro_mu << endl;
	cerr << "ALPHA = " << distro_alpha << endl; 

      // OUTPUT THE ORIGINAL HISTOGRAM
	cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
	for (size_t i = 0; i < counts_hist.size(); i++)
	  if (counts_hist[i] > 0)
	    cerr << i << '\t' << setprecision(16) << counts_hist[i] << endl;
      }

      vector<double> measure_moments;
  // mu_r = (r + 1)! n_{r+1} / n_1
      size_t indx = 1;
      while(counts_hist[indx] > 0  && indx <= counts_hist.size()){
	measure_moments.push_back(exp(gsl_sf_lnfact(indx)
				      + log(counts_hist[indx])
				      - log(counts_hist[1])));
	if(!std::isfinite(measure_moments.back())){
	  measure_moments.pop_back();
	  break;
	}
	indx++;
      }
      

      size_t n_points = std::min(num_points, static_cast<size_t>(floor(measure_moments.size()/2)));

      if(n_points != num_points && VERBOSE)
	cerr << "n_points = " << n_points << endl;      

      if(VERBOSE){
	cerr << "MOMENTS" << endl;
	for(size_t i = 0; i < measure_moments.size(); i++)
	  cerr << std::setprecision(16) << measure_moments[i] << endl;
	cerr << "OBSERVED_DISTINCT = " << distinct_reads << endl;
      }
    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////

      vector<double> quad_estimates;
      quadrature_bootstraps(VERBOSE, counts_hist, bootstraps, n_points,
			    tolerance, quad_estimates);
  
      double log_mean, log_lower_ci, log_upper_ci;
      log_mean_quad(VERBOSE, quad_estimates, ci_level,
		    log_mean, log_lower_ci, log_upper_ci);
    


      std::ofstream quad_of;
      if (!quad_outfile.empty()) quad_of.open(quad_outfile.c_str());
      std::ostream quad_out(quad_outfile.empty() ? std::cout.rdbuf() : quad_of.rdbuf());
      quad_out << "log_mean" << '\t' << "log_lower_" << ci_level << "ci" 
	       << '\t' << "log_upper_" << ci_level << "ci" << endl;
      quad_out << log_mean << '\t' << log_lower_ci << '\t'
	       << log_upper_ci << endl;
    
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
