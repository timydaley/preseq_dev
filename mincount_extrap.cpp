/*    lc_extrap:
 *
 *    Copyright (C) 2012 University of Southern California and
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

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <RNG.hpp>
#include <smithlab_os.hpp>

#include "continued_fraction.hpp"
#include "load_data_for_complexity.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;

using std::setw;
using std::fixed;
using std::setprecision;
using std::tr1::unordered_map;


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


static double
sample_count_reads_w_mincount(const gsl_rng *rng,
			      const vector<size_t> &full_umis,
			      const size_t sample_size,
			      const size_t mincount) {
  vector<size_t> sample_umis(sample_size);
  gsl_ran_choose(rng, (size_t *)&sample_umis.front(), sample_size,
		 (size_t *)&full_umis.front(), full_umis.size(), 
		 sizeof(size_t));

  double number_observed = 0.0;
  size_t current_count = 1;

  for (size_t i = 1; i < sample_umis.size(); i++){
    if(sample_umis[i] == sample_umis[i-1])
      current_count++;
    else{
      if(current_count >= mincount)
	number_observed++;
      current_count = 1;
    }
  }
  if(current_count >= mincount)
    number_observed++;

  return number_observed;
}


static bool
check_mincount_estimates_stability(const vector<double> &estimates,
				   const double max_change_per_time_step) {
  // make sure that the estimate is increasing in the time_step and
  // is below the initial distinct per step_size
  for (size_t i = 1; i < estimates.size(); ++i){
    if(!std::isfinite(estimates[i])){
      return false;
    }
    if ((estimates[i] < estimates[i - 1]) ||
	(estimates[i] - estimates[i - 1] > max_change_per_time_step)){
      return false;
    }
  }    

  return true;
}


void
estimates_bootstrap(const bool VERBOSE, const vector<double> &orig_hist, 
		    const size_t bootstraps, const size_t orig_max_terms, 
		    const int diagonal, const double step_size, 
		    const double max_extrapolation, const size_t mincount,
		    vector< vector<double> > &full_estimates) {
  // clear returning vectors
  full_estimates.clear();
  
  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand()); 


  const double initial_observed 
    = accumulate(orig_hist.begin() + mincount, orig_hist.end(), 0.0);


  vector<size_t> orig_hist_distinct_counts;
  vector<double> distinct_orig_hist;
  for (size_t i = 0; i < orig_hist.size(); i++){
    if (orig_hist[i] > 0) {
      orig_hist_distinct_counts.push_back(i);
      distinct_orig_hist.push_back(orig_hist[i]);
    }
  }
  
  double vals_sum = 0.0;
  for(size_t i = 0; i < orig_hist.size(); i++)
    vals_sum += i*orig_hist[i];
    
  for (size_t iter = 0; 
       (iter < 2*bootstraps && full_estimates.size() < bootstraps); ++iter) {
    
    vector<double> mincount_vector;
    vector<double> hist;
    resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);
    
    //resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();

    //construct umi vector to sample from
    vector<size_t> umis;
    size_t umi = 1;
    for(size_t i = 1; i < hist.size(); i++){
      for(size_t j = 0; j < hist[i]; j++){
	for(size_t k = 0; k < i; k++)
	  umis.push_back(umi);
	umi++;
      }
    }
    assert(umis.size() == static_cast<size_t>(vals_sum));

    // compute complexity curve by random sampling w/out replacement
    size_t upper_limit = static_cast<size_t>(vals_sum);
    size_t step = static_cast<size_t>(step_size);
    size_t sample = step;
    while(sample <= upper_limit){
      mincount_vector.push_back(sample_count_reads_w_mincount(rng, umis, sample, mincount));
      sample += step;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < hist.size() &&
	   hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;   

    size_t max_terms = std::min(orig_max_terms, counts_before_first_zero - mincount);
    max_terms = max_terms - (max_terms % 2 == 0);

    const ContinuedFractionApproximation 
      mincount_cfa(diagonal, max_terms);
    
    const ContinuedFraction 
      mincount_cf(mincount_cfa.optimal_cont_frac_mincount(hist, mincount, diagonal));


    //extrapolate the curve start
    if (mincount_cf.is_valid()){
      //  cerr << "valid" << endl;
      double sample_size = static_cast<double>(sample);
      while(sample_size <= max_extrapolation){
	double t = (sample_size - vals_sum)/vals_sum;
	assert(t >= 0.0);
	mincount_vector.push_back(initial_observed + t*mincount_cf(t));
	sample_size += step_size;
      }
      if(check_mincount_estimates_stability(mincount_vector, step_size)){
	full_estimates.push_back(mincount_vector);
	if (VERBOSE) cerr << '.';
      }
      else if(VERBOSE) cerr << '_';
    }
  else if(VERBOSE) cerr << '_';
  }
  if (VERBOSE)
    cerr << endl;
  if (full_estimates.size() < bootstraps)
    throw SMITHLABException("too many iterations, poor sample");
}


static void
return_median_and_ci(const vector<vector<double> > &estimates,
		     const double alpha, const double initial_distinct,
		     const double vals_sum, const double step_size, 
		     vector<double> &median_estimates,
		     vector<double> &lower_ci, vector<double> &upper_ci) {
  //  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
  assert(!estimates.empty());
  
  const size_t n_est = estimates.size();
  vector<double> estimates_row(estimates.size(), 0.0);
  for (size_t i = 0; i < estimates[0].size(); i++) {
    
    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = estimates[k][i];
    
    // sort to get confidence interval
    sort(estimates_row.begin(), estimates_row.end());
    const double curr_median = 
      gsl_stats_median_from_sorted_data(&estimates_row[0], 1, n_est);
    
    median_estimates.push_back(curr_median);
    // const double std_dev = sqrt(gsl_stats_variance(&estimates_row[0], 1, n_est));
    upper_ci.push_back(gsl_stats_quantile_from_sorted_data(&estimates_row[0], 1, n_est, 0.95));
    lower_ci.push_back(gsl_stats_quantile_from_sorted_data(&estimates_row[0], 1, n_est, 0.05));
  }
}

static bool
mincount_single_estimate(const bool VERBOSE, vector<double> &hist,
                       size_t orig_max_terms, const int diagonal,
                       const double step_size, 
                       const double max_extrapolation,
		       const size_t mincount,
                       vector<double> &mincount_estimate) {

  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand());


  mincount_estimate.clear();
  double vals_sum = 0.0;
  for(size_t i = 0; i < hist.size(); i++)
    vals_sum += i*hist[i];
  const double initial_observed 
    = accumulate(hist.begin() + mincount, hist.end(), 0.0);

  //construct umi vector to sample from
  vector<size_t> umis;
  size_t umi = 1;
  for(size_t i = 1; i < hist.size(); i++){
    for(size_t j = 0; j < hist[i]; j++){
      for(size_t k = 0; k < i; k++)
        umis.push_back(umi);
      umi++;
    }
  }
  assert(umis.size() == static_cast<size_t>(vals_sum));


    // compute complexity curve by random sampling w/out replacement
  size_t upper_limit = static_cast<size_t>(vals_sum);
  size_t step = static_cast<size_t>(step_size);
  size_t sample = step;
  while(sample <= upper_limit){
    mincount_estimate.push_back(sample_count_reads_w_mincount(rng, umis, sample, mincount));
    sample += step;
  }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
  size_t counts_before_first_zero = 1;
  while (counts_before_first_zero < hist.size() &&
	 hist[counts_before_first_zero] > 0)
    ++counts_before_first_zero;   

  size_t max_terms = std::min(orig_max_terms, counts_before_first_zero - mincount);
  max_terms = max_terms - (max_terms % 2 == 0);

  const ContinuedFractionApproximation 
    mincount_cfa(diagonal, max_terms);
    
  const ContinuedFraction 
    mincount_cf(mincount_cfa.optimal_cont_frac_mincount(hist, mincount, diagonal));


    //extrapolate the curve start
  if (mincount_cf.is_valid()){
      //  cerr << "valid" << endl;
    double sample_size = static_cast<double>(sample);
    while(sample_size <= max_extrapolation){
      double t = (sample_size - vals_sum)/vals_sum;
      assert(t >= 0.0);
      mincount_estimate.push_back(initial_observed + t*mincount_cf(t));
      sample_size += step_size;
    }
  }
  else{
    // FAIL!
    // lower_cf unacceptable, need to bootstrap to obtain estimates
    return false;
  }

  if (VERBOSE) {
    if(mincount_cf.offset_coeffs.size() > 0){
      cerr << "CF_OFFSET_COEFF_ESTIMATES" << endl;
      copy(mincount_cf.offset_coeffs.begin(), mincount_cf.offset_coeffs.end(),
           std::ostream_iterator<double>(cerr, "\n"));
    }
    if(mincount_cf.cf_coeffs.size() > 0){
      cerr << "CF_COEFF_ESTIMATES" << endl;
      copy(mincount_cf.cf_coeffs.begin(), mincount_cf.cf_coeffs.end(),
           std::ostream_iterator<double>(cerr, "\n"));
    }
  }

  // SUCCESS!!
  return true;
}


static void
write_predicted_curve(const string outfile, const double c_level,
		      const double step_size,
		      const vector<double> &median_mincount_estimates,
		      const vector<double> &mincount_lower_ci,
		      const vector<double> &mincount_upper_ci) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  out << "TOTAL_READS\tEXPECTED_DISTINCT\t"
      << "LOWER_" << 100*c_level << "%CI\t"
      << "UPPER_" << 100*c_level << "%CI" << endl;
  
  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);
  
  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < median_mincount_estimates.size(); ++i)
    out << (i + 1)*step_size << '\t' 
	<< median_mincount_estimates[i] << '\t'
	<< mincount_lower_ci[i] << '\t' << mincount_upper_ci[i] << endl;
}



int
main(const int argc, const char **argv) {

  try {
    
    //   const size_t MIN_REQUIRED_COUNTS = 8;

    /* FILES */
    string outfile;
    
    size_t orig_max_terms = 200;
    double max_extrapolation = 1.0e10;
    double step_size = 1e6;
    size_t bootstraps = 100;
    int diagonal = 0;
    double c_level = 0.95;
    size_t mincount = 2;
    
    /* FLAGS */
    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool SINGLE_ESTIMATE = false;

    
#ifdef HAVE_SAMTOOLS
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
#endif
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "mincount output file (default: stdout)",
		      false , outfile);
    opt_parse.add_opt("extrap",'e',"maximum extrapolation "
		      "(default: " + toa(max_extrapolation) + ")",
		      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size in extrapolations "
		      "(default: " + toa(step_size) + ")", 
		      false, step_size);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps "
		      "(default: " + toa(bootstraps) + "), ",
		      false, bootstraps);
    opt_parse.add_opt("cval", 'c', "level for confidence intervals "
		      "(default: " + toa(c_level) + ")", false, c_level);
    opt_parse.add_opt("mincount",'m',"minimum number of observations for a read to be counted "
		      "(default: " + toa(mincount) + "), ",
		      false, mincount); 
    //	opt_parse.add_opt("terms",'t',"maximum number of terms", false, 
    //	     orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false, VERBOSE);
#ifdef HAVE_SAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format", 
		      false, BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l', "maximum segment length when merging "
                      "paired end bam reads (default: "
                      + toa(MAX_SEGMENT_LENGTH) + ")",
                      false, MAX_SEGMENT_LENGTH);
#endif
    opt_parse.add_opt("pe", 'P', "input is paired end read file",
		      false, PAIRED_END);
    opt_parse.add_opt("vals", 'V', 
		      "input is a text file containing only the observed counts",
		      false, VALS_INPUT);
    opt_parse.add_opt("hist", 'H',
                      "input is a text file containing the observed histogram",
                      false, HIST_INPUT);
    opt_parse.add_opt("quick",'Q',
                      "quick mode, estimate mincount without bootstrapping for confidence intervals",
                      false, SINGLE_ESTIMATE);

    
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
    /******************************************************************/

        
    vector<double> counts_hist;
    size_t n_reads = 0;

    // LOAD VALUES
    if(HIST_INPUT){
      if(VERBOSE)
        cerr << "HIST_INPUT" << endl;
      n_reads = load_histogram(input_file_name, counts_hist);
    }
    else if(VALS_INPUT){
      if(VERBOSE)
        cerr << "VALS_INPUT" << endl;
      n_reads = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_SAMTOOLS
    else if (BAM_FORMAT_INPUT && PAIRED_END){
      if(VERBOSE)
        cerr << "PAIRED_END_BAM_INPUT" << endl;
      const size_t MAX_READS_TO_HOLD = 5000000;
      size_t n_paired = 0;
      size_t n_mates = 0;
      n_reads = load_counts_BAM_pe(VERBOSE, input_file_name, 
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
      n_reads = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if(PAIRED_END){
      if(VERBOSE)
        cerr << "PAIRED_END_BED_INPUT" << endl;
      n_reads = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else{ // default is single end bed file
      if(VERBOSE)
        cerr << "BED_INPUT" << endl;
      n_reads = load_counts_BED_se(input_file_name, counts_hist);
    }

    const size_t max_observed_count = counts_hist.size() - 1;
    //   const double distinct_reads = accumulate(counts_hist.begin(),
    //                                       counts_hist.end(), 0.0);

    // for large initial experiments need to adjust step size
    // otherwise small relative steps do not account for variance
    // in extrapolation
    if(step_size < (n_reads/20)){
      step_size = std::max(step_size, 
                           step_size*round(n_reads/(20*step_size)));
      if(VERBOSE)
        cerr << "ADJUSTED_STEP_SIZE = " << step_size << endl;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < counts_hist.size() &&
           counts_hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;

    orig_max_terms = std::min(orig_max_terms, counts_before_first_zero - 1);
    orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);

    const double initial_observed 
      = accumulate(counts_hist.begin() + mincount, counts_hist.end(), 0.0);

    double values_sum = 0.0;
    for(size_t i = 0; i < counts_hist.size(); i++)
	  values_sum += i*counts_hist[i];


    const size_t distinct_counts =
      static_cast<size_t>(std::count_if(counts_hist.begin(), counts_hist.end(),
                                        bind2nd(std::greater<double>(), 0.0)));
    if (VERBOSE)
      cerr << "TOTAL READS      = " << n_reads << endl
           << "INITIAL MINCOUNT = " << initial_observed << endl
           << "DISTINCT COUNTS  = " << distinct_counts << endl
           << "MAX COUNT        = " << max_observed_count << endl
           << "COUNTS OF 1      = " << counts_hist[1] << endl
           << "MAX TERMS        = " << orig_max_terms << endl;

    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
        if (counts_hist[i] > 0)
          cerr << i << '\t' << static_cast<size_t>(counts_hist[i]) << endl;
      cerr << endl;
    }


    if(SINGLE_ESTIMATE){
      vector<double> mincount_estimates;
      bool SINGLE_ESTIMATE_SUCCESS = 
	mincount_single_estimate(VERBOSE, counts_hist, 
				  orig_max_terms, diagonal, 
				  step_size, max_extrapolation, 
				  mincount,
				  mincount_estimates);

     // IF FAILURE, EXIT
      if(!SINGLE_ESTIMATE_SUCCESS)
        throw SMITHLABException("SINGLE ESTIMATE FAILED, NEED TO RUN "
                                "FULL MODE FOR ESTIMATES");

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS\tEXPECTED_DISTINCT" << endl;

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << 0 << '\t' << 0 << endl;
      for (size_t i = 0; i < mincount_estimates.size(); ++i)
        out << (i + 1)*step_size << '\t'
            << mincount_estimates[i] << endl;

    }
    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // BOOTSTRAPS

    else{
      if (VERBOSE) 
	cerr << "[BOOTSTRAP ESTIMATES]" << endl;
      
      vector<vector <double> > mincount_estimates;
      estimates_bootstrap(VERBOSE, counts_hist,  bootstraps, orig_max_terms,
			  diagonal, step_size, max_extrapolation, 
			  mincount, mincount_estimates);
      
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // BOOTSTRAPS
      if (VERBOSE)
	cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
      
      vector<double> median_mincount_estimates;
      vector<double> mincount_upper_ci, mincount_lower_ci;
      return_median_and_ci(mincount_estimates, 1.0 - c_level, 
			   initial_observed, values_sum, step_size,
			   median_mincount_estimates, 
			   mincount_lower_ci, mincount_upper_ci);

    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
      if (VERBOSE) 
	cerr << "[WRITING OUTPUT]" << endl;
    
      write_predicted_curve(outfile, c_level, step_size,
			    median_mincount_estimates,
			    mincount_lower_ci, mincount_upper_ci);
    }
      
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
