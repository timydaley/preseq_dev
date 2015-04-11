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
#include <queue>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <tr1/unordered_map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>


#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <smithlab_os.hpp>
#include <RNG.hpp>


#include "moment_sequence.hpp"
#include "load_data_for_complexity.hpp"


using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::fixed;
using std::setprecision;
using std::min;


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

static inline double
alpha_log_confint_multiplier(const double estimate,
                             const double variance, const double alpha) {
  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
  return exp(inv_norm_alpha*
             sqrt(log(1.0 + variance/pow(estimate, 2))));
}


static void
median_and_ci(const vector<double> &estimates,
              const double c_level,
              double &median_estimate,
              double &lower_ci_estimate,
              double &upper_ci_estimate){
  assert(!estimates.empty());
  const double alpha = 1.0 - c_level;
  const size_t n_est = estimates.size();
  vector<double> sorted_estimates(estimates);
  sort(sorted_estimates.begin(), sorted_estimates.end());
  median_estimate =
    gsl_stats_median_from_sorted_data(&sorted_estimates[0], 
                                      1, n_est);
  const double variance = 
    gsl_stats_variance(&sorted_estimates[0], 1, n_est);
  const double confint_mltr =
    alpha_log_confint_multiplier(median_estimate, variance, alpha);

  lower_ci_estimate = median_estimate/confint_mltr;
  upper_ci_estimate = median_estimate*confint_mltr;

}

void
log_mean(const bool VERBOSE,
	 const vector<double> &estimates,
	 const double c_level,
	 double &log_mean, 
	 double &log_lower_ci,
	 double &log_upper_ci){
  vector<double> log_estimates(estimates);
  for(size_t i = 0; i < log_estimates.size(); i++)
    log_estimates[i] = log(log_estimates[i]);

  log_mean = exp(gsl_stats_mean(&log_estimates[0], 1,
				log_estimates.size()) );

  double log_std_dev = std::sqrt(gsl_stats_variance(&log_estimates[0], 1, 
						    log_estimates.size()) );

  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv((1.0 - c_level)/2.0);
  log_lower_ci = exp(log(log_mean) - inv_norm_alpha*log_std_dev);
  log_upper_ci = exp(log(log_mean) + inv_norm_alpha*log_std_dev);
}

int
main(const int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    bool QUICK_MODE = false;

    string outfile;

#ifdef HAVE_SAMTOOLS
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
#endif

    size_t max_num_points = 5;
    double tolerance = 1e-20;
    size_t max_iter = 100;
    size_t bootstraps = 100;
    double c_level = 0.95;


    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[0]),
                           "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("max_num_points",'p',"maximum number of points in quadrature "
                      "estimates (default: " + toa(max_num_points) + ")",
                      false, max_num_points);
    opt_parse.add_opt("tolerance", 't', "numerical tolerance "
                      "(default: " + toa(tolerance) + ")",
		      false, tolerance);
    opt_parse.add_opt("bootstraps", 'n', "number of bootstraps "
                      "(default: " + toa(bootstraps) + ")",
		      false, bootstraps);
    opt_parse.add_opt("clevel", 'c', "level for confidence intervals "
                      "(default: " + toa(c_level) + ")", false, c_level);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    opt_parse.add_opt("pe", 'P', "input is paired end read file",
                      false, PAIRED_END);
    opt_parse.add_opt("hist", 'H',
                      "input is a text file containing the observed histogram",
                      false, HIST_INPUT);
    opt_parse.add_opt("vals", 'V',
                      "input is a text file containing only the observed counts",
                      false, VALS_INPUT);
#ifdef HAVE_SAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format",
                      false, BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l', "maximum segment length when merging "
                      "paired end bam reads (default: "
                      + toa(MAX_SEGMENT_LENGTH) + ")",
                      false, MAX_SEGMENT_LENGTH);
#endif
    opt_parse.add_opt("quick", 'q', "quick mode, estimates without bootstrapping",
		      false, QUICK_MODE);

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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file_name = leftover_args.front();
    // ****************************************************************
 
   vector<double> counts_hist;
    size_t n_obs = 0;

    // LOAD VALUES
    if(HIST_INPUT){
      if(VERBOSE)
        cerr << "HIST_INPUT" << endl;
      n_obs = load_histogram(input_file_name, counts_hist);
    }
    else if(VALS_INPUT){
      if(VERBOSE)
        cerr << "VALS_INPUT" << endl;
      n_obs = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_SAMTOOLS
    else if (BAM_FORMAT_INPUT && PAIRED_END){
      if(VERBOSE)
        cerr << "PAIRED_END_BAM_INPUT" << endl;
      const size_t MAX_READS_TO_HOLD = 5000000;
      size_t n_paired = 0;
      size_t n_mates = 0;
      n_obs = load_counts_BAM_pe(VERBOSE, input_file_name, 
                                   MAX_SEGMENT_LENGTH, 
                                   MAX_READS_TO_HOLD, n_paired, 
                                   n_mates, counts_hist);
      if(VERBOSE){
        cerr << "MERGED PAIRED END READS = " << n_paired << endl;
        cerr << "MATES PROCESSED = " << n_mates << endl;
      }
    }
    else if(BAM_FORMAT_INPUT){
      if(VERBOSE)
        cerr << "BAM_INPUT" << endl;
      n_obs = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if(PAIRED_END){
      if(VERBOSE)
        cerr << "PAIRED_END_BED_INPUT" << endl;
      n_obs = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else{ // default is single end bed file
      if(VERBOSE)
        cerr << "BED_INPUT" << endl;
      n_obs = load_counts_BED_se(input_file_name, counts_hist);
    }

    const double distinct_obs = accumulate(counts_hist.begin(), 
					   counts_hist.end(), 0.0);


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
    
    
    if (VERBOSE){
      cerr << "TOTAL OBSERVATIONS     = " << n_obs << endl
           << "DISTINCT OBSERVATIONS  = " << distinct_obs << endl
           << "MAX COUNT              = " << counts_hist.size() - 1 << endl;

      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << setprecision(16) << counts_hist[i] << endl;

      cerr << "OBSERVED MOMENTS" << endl;
      for(size_t i = 0; i < min(measure_moments.size(),
				2*max_num_points); i++)
	cerr << std::setprecision(16) << measure_moments[i] << endl;  
    }
    
    ////////////////////////////////////////////////////////////////////
    // calculate lower bound
    

    

    if(QUICK_MODE){
      if(measure_moments.size() < 2*max_num_points)
	max_num_points = static_cast<size_t>(floor(measure_moments.size()/2));
      else
	measure_moments.resize(2*max_num_points);
      size_t n_points = 0;
      n_points = ensure_pos_def_mom_seq(measure_moments, tolerance, VERBOSE);
      if(VERBOSE)
	cerr << "n_points = " << n_points << endl;    

      MomentSequence obs_mom_seq(measure_moments);
    
      if(VERBOSE){
	cerr << "alpha = ";
	for(size_t k = 0; k < obs_mom_seq.alpha.size(); k++)
	  cerr << obs_mom_seq.alpha[k] << ", ";
	cerr << endl;

	cerr << "beta = ";
	for(size_t k = 0; k < obs_mom_seq.beta.size(); k++)
	  cerr << obs_mom_seq.beta[k] << ", ";
	cerr << endl;
      }
    
      vector<double> points, weights;

      obs_mom_seq.Lower_quadrature_rules(VERBOSE, n_points, tolerance,
					 max_iter, points, weights);

      const double weights_sum = accumulate(weights.begin(), weights.end(), 0.0);
      if(weights_sum != 1.0){
	if(VERBOSE)
	  cerr << "weights sum = " << weights_sum << endl;
	for(size_t i = 0; i < weights.size(); i++)
	  weights[i] = weights[i]/weights_sum;
      }

      if(VERBOSE){
	cerr << "points = ";
	for(size_t i = 0; i < points.size(); i++)
	  cerr << setprecision(16) << points[i] << ", ";
	cerr << endl;

	cerr << "weights = ";
	for(size_t i = 0; i < weights.size(); i++)
	  cerr << setprecision(16) << weights[i] << ", ";
	cerr << endl;
      }
    
      double estimated_unobs = 0.0;
    
      for(size_t i = 0; i < weights.size(); i++)
	estimated_unobs += counts_hist[1]*weights[i]/points[i];

      if(estimated_unobs > 0.0)
	estimated_unobs += distinct_obs;
      else{
	estimated_unobs = distinct_obs;
	n_points = 0;
      }
    
    
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "quadrature_estimated_unobs" << '\t' << "n_points" << endl;
      out << estimated_unobs << '\t' << n_points << endl;
    
    }
    else{
      vector<double> quad_estimates;

  //setup rng
      srand(time(0) + getpid());
      gsl_rng_env_setup();
      gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
      gsl_rng_set(rng, rand());

      // hist may be sparse, to speed up bootstrapping
      // sample only from positive entries
      vector<size_t> counts_hist_distinct_counts;
      vector<double> distinct_counts_hist;
      for (size_t i = 0; i < counts_hist.size(); i++){
	if (counts_hist[i] > 0) {
	  counts_hist_distinct_counts.push_back(i);
	  distinct_counts_hist.push_back(counts_hist[i]);
	}
      }

      for(size_t iter = 0; 
	  iter < max_iter && quad_estimates.size() < bootstraps; 
	  ++iter){

	vector<double> sample_hist;
	resample_hist(rng, counts_hist_distinct_counts, 
		      distinct_counts_hist, sample_hist);

	const double sampled_distinct = accumulate(sample_hist.begin(), sample_hist.end(), 0.0);
	// initialize moments, 0th moment is 1
	vector<double> bootstrap_moments(1, 1.0);
	// moments[r] = (r + 1)! n_{r+1} / n_1
	for(size_t i = 0; i < 2*max_num_points + 1; i++)
	  bootstrap_moments.push_back(exp(gsl_sf_lnfact(i + 2) 
					  + log(sample_hist[i + 2])
					  - log(sample_hist[1])) );

	size_t n_points = 0;
	n_points = ensure_pos_def_mom_seq(bootstrap_moments, tolerance, VERBOSE);


	MomentSequence bootstrap_mom_seq(bootstrap_moments);

   	vector<double> points, weights;
	bootstrap_mom_seq.Lower_quadrature_rules(VERBOSE, n_points, tolerance,
						 max_iter, points, weights);

	const double weights_sum = accumulate(weights.begin(), weights.end(), 0.0);
	if(weights_sum != 1.0){
	  for(size_t i = 0; i < weights.size(); i++)
	    weights[i] = weights[i]/weights_sum;
	}

	double estimated_unobs = 0.0;
    
	for(size_t i = 0; i < weights.size(); i++)
	  estimated_unobs += counts_hist[1]*weights[i]/points[i];

	if(estimated_unobs > 0.0)
	  estimated_unobs += sampled_distinct;
	else{
	  estimated_unobs = sampled_distinct;
	  n_points = 0;
	}

	quad_estimates.push_back(estimated_unobs);
      }

      double median_estimate, log_mean_estimate, lower_log_ci, upper_log_ci;

      log_mean(VERBOSE, quad_estimates, c_level, log_mean_estimate, 
	       lower_log_ci, upper_log_ci);
      median_and_ci(quad_estimates, c_level, median_estimate,
		    lower_log_ci, upper_log_ci);

     std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "median_estimated_unobs" << '\t' 
	  << "log_mean_estimated_unobs" << '\t'
	  << "log_lower_ci" << '\t'
	  << "log_upper_ci" << endl;
      out << median_estimate << '\t' 
	  << log_mean_estimate << '\t'
	  << lower_log_ci << '\t'
	  << upper_log_ci << endl;


    }
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // done trying

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
