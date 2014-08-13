/*    saturation_extrap:
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith and Timothy Daley
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


#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>
#include <queue>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <tr1/unordered_map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>


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
using std::min;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::ifstream;
using std::isfinite;
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

// sample umis without replacement to interpolate complexity
static double
sample_count_singletons(const gsl_rng *rng,
			const vector<size_t> &full_umis,
			const size_t sample_size) {
    vector<size_t> sample_umis(sample_size);
    gsl_ran_choose(rng, (size_t *)&sample_umis.front(), sample_size,
                   (size_t *)&full_umis.front(), full_umis.size(),
                   sizeof(size_t));
    double n_singletons = 0.0;
    if(sample_umis[1] != sample_umis[0])
      ++n_singletons;
    for (size_t i = 2; i < sample_umis.size(); i++)
      if((sample_umis[i] != sample_umis[i-1])
	 && (sample_umis[i-2] != sample_umis[i-1]))
            ++n_singletons;
    
    return n_singletons;
}


static inline bool
check_saturation_estimates(const vector<double> estimates){
  if(estimates.empty())
    return false;

  // make sure estimates are decreasing and
  // between 0 & 1
  if(estimates[0] >= 1.0 || estimates[0] < 0.0)
    return false;

  for(size_t i = 1; i < estimates.size(); i++)
    if(estimates[i] > estimates[i-1] ||
       estimates[i] >= 1.0 ||
       estimates[i] < 0.0) 
      return false;
  
  return true;
}

void
bootstrap_saturation_deriv(const bool VERBOSE, const vector<double> &orig_hist,
			   const size_t bootstraps, const size_t orig_max_terms, 
			   const int diagonal, const double step_size, 
			   const double max_extrapolation, 
			   vector< vector<double> > &full_estimates){
  full_estimates.clear();

    //setup rng
    srand(time(0) + getpid());
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, rand());
    
    double vals_sum = 0.0;
    for(size_t i = 0; i < orig_hist.size(); i++)
        vals_sum += orig_hist[i]*i;    
    
    vector<size_t> orig_hist_distinct_counts;
    vector<double> distinct_orig_hist;
    for(size_t i = 0; i < orig_hist.size(); i++){
        if(orig_hist[i] > 0){
            orig_hist_distinct_counts.push_back(i);
            distinct_orig_hist.push_back(orig_hist[i]);
        }
    }

  for (size_t iter = 0; 
       iter < 10*bootstraps && full_estimates.size() < bootstraps; ++iter) {

    vector<double> yield_vector;
    vector<double> hist;
    resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);
        
    double sample_vals_sum = 0.0;
    for(size_t i = 0; i < hist.size(); i++)
      sample_vals_sum += i*hist[i];
        
    const double sample_max_val = max_extrapolation/sample_vals_sum;
    const double sample_val_step = step_size/sample_vals_sum;
        
        //resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();
    

      // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < hist.size() &&
	   hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;
        
    size_t max_terms = std::min(orig_max_terms, counts_before_first_zero - 1);
        // refit curve for lower bound (degree of approx is 1 less than
        // max_terms)
    max_terms = max_terms - (max_terms % 2 == 1);
        
        //refit curve for lower bound
    const ContinuedFractionApproximation
      lower_cfa(diagonal, max_terms);
        
    const ContinuedFraction
      lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));
        
    vector<double> saturation_estimates;
    if(lower_cf.is_valid()){

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
        assert(umis.size() == static_cast<size_t>(sample_vals_sum));
        
        // interpolate complexity curve by random sampling w/out replacement
        size_t upper_limit = static_cast<size_t>(sample_vals_sum);
        size_t step = static_cast<size_t>(step_size);
        size_t sample = step;
        while(sample < upper_limit){
            saturation_estimates.push_back(sample_count_singletons(rng, umis, sample)/sample);
            sample += step;
        }

	lower_cf.extrapolate_yield_deriv(hist, vals_sum, 
					 static_cast<double>(sample)/sample_vals_sum,
					 sample_max_val, sample_val_step, 
					 saturation_estimates);

    }
    else if(VERBOSE)
      cerr << "not_valid" << endl;

    if(check_saturation_estimates(saturation_estimates)){
      full_estimates.push_back(saturation_estimates);
      if(VERBOSE)
	cerr << ".";
    }
    else if(VERBOSE){
      cerr << "_" << endl;
    }

  }
}

static bool
extrap_single_estimate(const bool VERBOSE, vector<double> &hist,
                       size_t max_terms, const int diagonal,
                       const double step_size, 
                       const double max_extrapolation,
                       vector<double> &saturation_estimates) {

  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand());


  saturation_estimates.clear();
  double vals_sum = 0.0;
  for(size_t i = 0; i < hist.size(); i++)
    vals_sum += i*hist[i];

  const double max_val = max_extrapolation/vals_sum;
  const double val_step = step_size/vals_sum;

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

  // interpolate complexity curve by random sampling w/out replacement
  size_t upper_limit = static_cast<size_t>(vals_sum);
  size_t step = static_cast<size_t>(step_size);
  size_t sample = step;
  while (sample < upper_limit){
    saturation_estimates.push_back(sample_count_singletons(rng, umis, sample)/sample);
    sample += step;
  }

  // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
  size_t counts_before_first_zero = 1;
  while (counts_before_first_zero < hist.size() &&
         hist[counts_before_first_zero] > 0)
    ++counts_before_first_zero;


  // Ensure we are not using a zero term
  max_terms = std::min(max_terms, counts_before_first_zero - 1);

  // refit curve for lower bound (degree of approx is 1 less than
  // max_terms)
  max_terms = max_terms - (max_terms % 2 == 1);

  const ContinuedFractionApproximation
    lower_cfa(diagonal, max_terms);

  const ContinuedFraction
    lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));

  // extrapolate curve
  if (lower_cf.is_valid()){
	lower_cf.extrapolate_yield_deriv(hist, vals_sum, 
					 static_cast<double>(sample)/vals_sum,
					 max_val, val_step, 
					 saturation_estimates);
  }
  else{
    // FAIL!
    // lower_cf unacceptable, need to bootstrap to obtain estimates
    return false;
  }

  if (VERBOSE) {
    if(lower_cf.offset_coeffs.size() > 0){
      cerr << "CF_OFFSET_COEFF_ESTIMATES" << endl;
      copy(lower_cf.offset_coeffs.begin(), lower_cf.offset_coeffs.end(),
           std::ostream_iterator<double>(cerr, "\n"));
    }
    if(lower_cf.cf_coeffs.size() > 0){
      cerr << "CF_COEFF_ESTIMATES" << endl;
      copy(lower_cf.cf_coeffs.begin(), lower_cf.cf_coeffs.end(),
           std::ostream_iterator<double>(cerr, "\n"));
    }
  }

  // SUCCESS!!
  return true;
}


static double
compute_var(const vector<double> &estimates,
	    const double mean) {
  double variance = 0.0;
  for(size_t i = 0; i < estimates.size(); i++)
    variance += (estimates[i] - mean)*(estimates[i] - mean)/estimates.size();
  return variance; 
}


static void
compute_mean_and_alphaCI(const vector< vector<double> > &full_estimates,
			 const double alpha, 
			 vector<double> &mean_estimates,
			 vector<double> &lower_alphaCI,
			 vector<double> &upper_alphaCI){
  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);

  for(size_t i = 0; i < full_estimates[0].size(); i++) {
    // estimates is in wrong order, work locally on const val
    vector<double> log_estimates_row(full_estimates.size(), 0.0);

    for(size_t k = 0; k < log_estimates_row.size(); ++k)
      log_estimates_row[k] = log(full_estimates[k][i]);

    mean_estimates.push_back(exp(accumulate(log_estimates_row.begin(),
					    log_estimates_row.end(), 0.0)/
				 log_estimates_row.size()));
    const double variance = compute_var(log_estimates_row,
					log(mean_estimates.back()));

    // log confidence intervals
    upper_alphaCI.push_back(exp(log(mean_estimates.back())
				    + inv_norm_alpha*sqrt(variance)));
    lower_alphaCI.push_back(exp(log(mean_estimates.back())
				    - inv_norm_alpha*sqrt(variance)));

  }

}

static double
GoodToulmin2xExtrap(const vector<double> &counts_hist){
    double two_fold_extrap = 0.0;
    for(size_t i = 0; i < counts_hist.size(); i++)
        two_fold_extrap += pow(-1.0, i + 1)*counts_hist[i];
    
    return two_fold_extrap;
}


int
main(const int argc, const char **argv) {

  try {
      const size_t MIN_REQUIRED_COUNTS = 4;
        
        /* FILES */
        string outfile;
        
        size_t orig_max_terms = 100;
        double max_extrapolation = 1.0e10;
        
        // AS: this step size issue needs to be addressed
        double step_size = 1e6;
        // double read_step_size = 1e7;
        
        size_t MAX_SEGMENT_LENGTH = 10000;
        //size_t max_width = 1000;
        size_t bootstraps = 100;
        int diagonal = -1;
        //size_t bin_size = 20;
        double c_level = 0.95;
        //double tolerance = 1e-20;
        //size_t max_iter = 0;
        double dupl_level = 0.5;
        //double reads_per_base = 2.0;
        //double fixed_fold = 20;
        
        /* FLAGS */
        //size_t MODE = 0;
        //bool NO_SEQUENCE = false;
        bool VERBOSE = false;
        bool VALS_INPUT = false;
        bool PAIRED_END = false;
        bool HIST_INPUT = false;
        bool SINGLE_ESTIMATE = false;
        
#ifdef HAVE_SAMTOOLS
        bool BAM_FORMAT_INPUT = false;
#endif
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
	OptionParser opt_parse(strip_path(argv[0]),
			       "", "<sorted-bed-file>");
	opt_parse.add_opt("output", 'o', "saturation output file (default: stdout)",
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
	opt_parse.add_opt("dupl_level",'d', "fraction of duplicate to predict "
			  "(default: " + toa(dupl_level) + ")",
			  false, dupl_level);
	opt_parse.add_opt("terms",'x',"maximum number of terms", false,
			  orig_max_terms);
            //    opt_parse.add_opt("tol",'t', "numerical tolerance", false, tolerance);
            //    opt_parse.add_opt("max_iter",'i', "maximum number of iteration",
            //		      false, max_iter);
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
			  "quick mode, estimate saturation without bootstrapping for confidence intervals",
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
    const double distinct_reads = accumulate(counts_hist.begin(),
                                             counts_hist.end(), 0.0);

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


    const size_t distinct_counts =
      static_cast<size_t>(std::count_if(counts_hist.begin(), counts_hist.end(),
                                        bind2nd(std::greater<double>(), 0.0)));
    if (VERBOSE)
      cerr << "TOTAL READS     = " << n_reads << endl
           << "DISTINCT READS  = " << distinct_reads << endl
           << "DISTINCT COUNTS = " << distinct_counts << endl
           << "MAX COUNT       = " << max_observed_count << endl
           << "COUNTS OF 1     = " << counts_hist[1] << endl
           << "MAX TERMS       = " << orig_max_terms << endl;

    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
        if (counts_hist[i] > 0)
          cerr << i << '\t' << static_cast<size_t>(counts_hist[i]) << endl;
      cerr << endl;
    }

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(counts_hist);
    if(two_fold_extrap < 0.0)
      throw SMITHLABException("Library expected to saturate in doubling of "
                              "size, unable to extrapolate");


    size_t total_reads = 0;
    for(size_t i = 0; i < counts_hist.size(); i++){
      total_reads += i*counts_hist[i];
    }
    //assert(total_reads == n_reads);

    // catch if all reads are distinct
    if (orig_max_terms < MIN_REQUIRED_COUNTS)
      throw SMITHLABException("max count before zero is les than min required "
                              "count (4), sample not sufficiently deep or "
                              "duplicates removed");


    if(SINGLE_ESTIMATE){
      vector<double> saturation_estimates;
      extrap_single_estimate(VERBOSE, counts_hist, orig_max_terms, diagonal,
			     step_size, max_extrapolation, saturation_estimates); 


      if(VERBOSE)
	cerr << "[WRITING OUTPUT]" << endl;

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS" << '\t' 
	  << "SATURATION" << endl;
    
      out << 0 << '\t' << 1.0 << endl;
      for (size_t i = 0; i < saturation_estimates.size(); ++i)
	out << (i + 1)*step_size << '\t' 
	    << saturation_estimates[i] << endl;

    
    }
    else{
      if(VERBOSE)
	cerr << "COMPUTING SATURATION" << endl;

      vector< vector<double> > full_deriv_estimates;
      bootstrap_saturation_deriv(VERBOSE, counts_hist,  
			   bootstraps, orig_max_terms, diagonal, 
			   step_size, max_extrapolation, 
			   full_deriv_estimates);



      if(VERBOSE)
	cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;

      vector<double> mean_deriv_estimates, upper_alphaCI_deriv, lower_alphaCI_deriv;
      compute_mean_and_alphaCI(full_deriv_estimates, 1.0 - c_level,
			       mean_deriv_estimates, lower_alphaCI_deriv,
			       upper_alphaCI_deriv);


      if(VERBOSE)
	cerr << "[WRITING OUTPUT]" << endl;

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS" << '\t' 
	  << "SATURATION" << '\t' 
	  << "LOWER_" << 100*c_level << "%CI" << '\t' 
	  << "UPPER_" << 100*c_level << "%CI" << endl;
    
      out << 0 << '\t' << 1.0 << '\t' << 1.0 << '\t' << 1.0 << endl;
      for (size_t i = 0; i < mean_deriv_estimates.size(); ++i)
	out << (i + 1)*step_size << '\t' 
	    << mean_deriv_estimates[i] << '\t'
	    << lower_alphaCI_deriv[i] << '\t' 
	    << upper_alphaCI_deriv[i] << endl;

    
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
