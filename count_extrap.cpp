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

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;

using std::setw;
using std::fixed;
using std::setprecision;
using std::tr1::unordered_map;

/*
 * This code is used to deal with read data in BAM format.
 */
#ifdef HAVE_BAMTOOLS
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;

static SimpleGenomicRegion
BamToSimpleGenomicRegion(const unordered_map<size_t, string> &chrom_lookup,
			 const BamAlignment &ba) {
  const unordered_map<size_t, string>::const_iterator
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));
  
  const string chrom = the_chrom->second;
  const size_t start = ba.Position;
  const size_t end = start + ba.Length;
  
  return SimpleGenomicRegion(chrom, start, end);
}


static size_t
load_values_BAM(const string &input_file_name, vector<double> &values) {
  
  BamReader reader;
  reader.Open(input_file_name);

  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;

  size_t n_reads = 1;
  values.push_back(1.0);

  SimpleGenomicRegion prev;
  BamAlignment bam;
  while (reader.GetNextAlignment(bam)) {
    SimpleGenomicRegion r(BamToSimpleGenomicRegion(chrom_lookup, bam));
    if (r.same_chrom(prev) && r.get_start() < prev.get_start())
      throw SMITHLABException("locations unsorted in: " + input_file_name);
    
    if (!r.same_chrom(prev) || r.get_start() != prev.get_start())
      values.push_back(1.0);
    else values.back()++;
    ++n_reads;
    prev.swap(r);
  }
  reader.Close();
  return n_reads;
}
#endif


static size_t
load_values(const string input_file_name, vector<double> &values) {
  
  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw SMITHLABException("problem opening file: " + input_file_name);
  
  SimpleGenomicRegion r, prev;
  if (!(in >> prev))
    throw SMITHLABException("problem reading from: " + input_file_name);
  
  size_t n_reads = 1;
  values.push_back(1.0);
  while (in >> r) {
    if (r < prev)
      throw SMITHLABException("locations unsorted in: " + input_file_name);
    if (!r.same_chrom(prev) || r.get_start() != prev.get_start())
      values.push_back(1.0);
    else values.back()++;
    ++n_reads;
    prev.swap(r);
  }
  return n_reads;
}


void
resample_hist(const gsl_rng *rng, const vector<double> &vals_hist,
	      const double total_sampled_reads,
	      double expected_sample_size,
	      vector<double> &sample_hist) {
  
  const size_t hist_size = vals_hist.size();
  const double vals_mean = total_sampled_reads/expected_sample_size;
  
  sample_hist = vector<double>(hist_size, 0.0);
  vector<unsigned int> curr_sample(hist_size);
  double remaining = total_sampled_reads;
  
  while (remaining > 0) {
    
    // get a new sample
    expected_sample_size = max(1.0, (remaining/vals_mean)/2.0);
    gsl_ran_multinomial(rng, hist_size, 
			static_cast<unsigned int>(expected_sample_size),
			&vals_hist.front(), &curr_sample.front());
    
    // see how much we got
    double inc = 0.0;
    for (size_t i = 0; i < hist_size; ++i)
      inc += i*curr_sample[i];
    
    // only add to histogram if sampled reads < remaining reads
    if (inc <= remaining) {
      for (size_t i = 0; i < hist_size; i++)
	sample_hist[i] += static_cast<double>(curr_sample[i]);
      // update the amount we still need to get
      remaining -= inc;
    }
  }
}


static bool
check_count_estimates(const vector<double> &estimates,
		      const size_t count) {
  if(estimates.size() == 0)
    return false;

  size_t number_modes = 0;
  for(size_t i = 1; i < estimates.size(); ++i){
    // make sure estimates are in bounds
    if(estimates[i] < 0.0 || estimates[i] > 3.2e9/count || !finite(estimates[i]))
      return false;
    // count modes by detecting change in sign of derivative
    if(i < estimates.size() - 1 &&
       ((estimates[i] - estimates[i - 1])
	*(estimates[i + 1] - estimates[i]) < 0))
      number_modes++;
  }
  // check unimodality
  if(number_modes > 1)
    return false;
 
  return true;
}

void
estimates_bootstrap(const bool VERBOSE, const vector<double> &orig_values, 
		    const size_t bootstraps, const size_t orig_max_terms, 
		    const int diagonal, const double step_size, 
		    const double max_extrapolation, 
		    const double max_val, const double val_step,
		    const size_t count,
		    vector< vector<double> > &count_estimates) {
  // clear returning vectors
  count_estimates.clear();
  
  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand()); 

  const size_t max_observed_count = 
    static_cast<size_t>(*std::max_element(orig_values.begin(), 
					  orig_values.end()));
    
  vector<double> orig_hist(max_observed_count + 1, 0.0);
  for (size_t i = 0; i < orig_values.size(); ++i)
    ++orig_hist[static_cast<size_t>(orig_values[i])];
  
  const double vals_sum = accumulate(orig_values.begin(), orig_values.end(), 0.0);
  
  for (size_t iter = 0; 
       (iter < (count + 1)*bootstraps && count_estimates.size() < bootstraps); 
       ++iter) {
    
    vector<double> hist;
    resample_hist(rng, orig_hist, vals_sum,
		  static_cast<double>(orig_values.size()), hist);
    
    //resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();
    
    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < hist.size() &&
	   hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;
    
    size_t max_terms = std::min(orig_max_terms, counts_before_first_zero - count);
    // refit curve for lower bound (degree of approx is 1 less than
    // max_terms)
    max_terms = max_terms - (max_terms % 2 == 0);
    
    //refit curve for lower bound
    const ContinuedFractionApproximation 
      lower_cfa(diagonal, max_terms, step_size, max_extrapolation);
    

    const ContinuedFraction 
      lower_cf(lower_cfa.optimal_cont_frac_count(hist, count));

    vector<double> count_vector;
    vector<double> sat_vector;
    if (lower_cf.is_valid()) 
      lower_cf.extrapolate_count(hist, max_val, val_step, count, count_vector);
    
    // SANITY CHECK    
    if (check_count_estimates(count_vector, count)) {
      count_estimates.push_back(count_vector);
      if (VERBOSE) cerr << '.';
    }
    else if (VERBOSE) 
      cerr << '_';
    
  }
  if (VERBOSE)
    cerr << endl;
  if (count_estimates.size() < bootstraps)
    throw SMITHLABException("too many iterations, poor sample");
}

static inline double
alpha_log_confint_multiplier(const double estimate,
			     const double initial_count,
			     const double variance, const double alpha) {
  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
  return exp(inv_norm_alpha*
	     sqrt(log(1.0 + variance/pow(estimate - initial_count, 2))));
}


static void
return_median_and_ci(const vector<vector<double> > &estimates,
		     const double alpha, const double initial_count,
		     vector<double> &median_estimates,
		     vector<double> &lower_ci, vector<double> &upper_ci) {
  
  assert(!estimates.empty());
  
  const size_t n_est = estimates.size();
  vector<double> estimates_row(estimates.size(), 0.0);
  median_estimates.push_back(initial_count);
  upper_ci.push_back(initial_count);
  lower_ci.push_back(initial_count);
  for (size_t i = 1; i < estimates[0].size(); i++) {
    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = estimates[k][i];
    
    // sort to get confidence interval
    sort(estimates_row.begin(), estimates_row.end());
    const double curr_median = 
      gsl_stats_median_from_sorted_data(&estimates_row[0], 1, n_est);
    
    median_estimates.push_back(curr_median);
    const double variance = gsl_stats_variance(&estimates_row[0], 1, n_est);
    cerr << "variance = " << variance << endl;
    if(!finite(variance)){
      cerr << "variance not finite" << endl;
      for(size_t j = 0; j < n_est; j++)
	cerr << estimates[j][i] << "\t";
      cerr << endl;
    }
    assert(finite(variance));
    if(!finite(curr_median)){
      cerr << "estimate not finite" << endl;
      for(size_t j = 0; j < n_est; j++)
	cerr << estimates[j][i] << "\t";
      cerr << endl;
    }

    const double confint_mltr = 
      alpha_log_confint_multiplier(curr_median, initial_count, 
				   variance, alpha);

    upper_ci.push_back(initial_count + 
		       (curr_median - initial_count)*confint_mltr);
    lower_ci.push_back(initial_count + 
		       (curr_median - initial_count)/confint_mltr);
  }
}


static void
write_predicted_curve(const string outfile, const double values_sum,
		      const double c_level, const double val_step,
		      const vector<double> &median_count_estimates,
		      const vector<double> &count_lower_ci,
		      const vector<double> &count_upper_ci) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  out << "TOTAL_READS\tEXPECTED_COUNT\t"
      << "LOWER_" << 100*c_level << "%CI\t"
      << "UPPER_" << 100*c_level << "%CI" << endl;
  
  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);
  
  double val = 0.0;
  for (size_t i = 0; i < median_count_estimates.size(); ++i, val += val_step)
    out << (val + 1.0)*values_sum << '\t' 
	<< median_count_estimates[i] << '\t'
	<< count_lower_ci[i] << '\t' << count_upper_ci[i] << endl;
}

/*
static void
write_estimates_outfile(const string outfile,
			const vector< vector<double> > &estimates,
			const double val_step, const double values_sum){

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  double val = 0.0;
  for(size_t i = 0; i < estimates[0].size(); i++, val += val_step){
    out << (val + 1.0)*values_sum << '\t';
    for(size_t j = 0; j < estimates.size(); j++)
      out << estimates[j][i] << "\t";
    out << endl;
  }
}
*/

int
main(const int argc, const char **argv) {

  try {
    
    const size_t MIN_REQUIRED_COUNTS = 8;

    /* FILES */
    string outfile;
    
    size_t orig_max_terms = 100;
    double max_extrapolation = 1.0e10;
    double step_size = 1e6;
    size_t bootstraps = 100;
    double c_level = 0.95;
    size_t count = 1;
    
    /* FLAGS */
    bool VERBOSE = false;
    // bool SMOOTH_HISTOGRAM = false;	
    
#ifdef HAVE_BAMTOOLS
    bool BAM_FORMAT_INPUT = false;
#endif
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "count output file (default: stdout)",
		      false , outfile);
    opt_parse.add_opt("count",'l', 
		      "count to extrapolate "
		      "(default:1)",
		      false, count);
    opt_parse.add_opt("extrap",'e',"maximum extrapolation "
		      "(default: " + toa(max_extrapolation) + ")",
		      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size in extrapolations "
		      "(default: " + toa(step_size) + ")", 
		      false, step_size);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps "
		      "(default: " + toa(bootstraps) + "), ",
		      false, bootstraps);
    opt_parse.add_opt("cval", 'i', "level for confidence intervals "
		      "(default: " + toa(c_level) + ")", false, c_level);
    //	opt_parse.add_opt("terms",'t',"maximum number of terms", false, 
    //	     orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false, VERBOSE);

#ifdef HAVE_BAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format", 
		      false, BAM_FORMAT_INPUT);
#endif
    
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
    const int diagonal = -count;

    vector<double> values;
#ifdef HAVE_BAMTOOLS
    if (BAM_FORMAT_INPUT)
      load_values_BAM(input_file_name, values);
    else
#endif
      load_values(input_file_name, values);
    
    // JUST A SANITY CHECK
    const size_t values_sum = 
      static_cast<size_t>(accumulate(values.begin(), values.end(), 0.0));
    
    const double max_val = max_extrapolation/static_cast<double>(values_sum);
    const double val_step = step_size/static_cast<double>(values_sum);
    
    const size_t max_observed_count = 
      static_cast<size_t>(*std::max_element(values.begin(), values.end()));

    // catch if all reads are distinct
    if (max_observed_count < MIN_REQUIRED_COUNTS)
      throw SMITHLABException("sample appears too uniform");
    
    // BUILD THE HISTOGRAM
    vector<double> counts_hist(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_hist[static_cast<size_t>(values[i])];

    const double initial_count = counts_hist[count];
    
    const size_t distinct_counts = 
      static_cast<size_t>(std::count_if(counts_hist.begin(), counts_hist.end(),
					bind2nd(std::greater<double>(), 0.0)));
    if (VERBOSE)
      cerr << "TOTAL READS       = " << values_sum << endl
	   << "DISTINCT READS    = " << values.size() << endl
	   << "INITIAL COUNT OF " << count << "= " << initial_count << endl
	   << "DISTINCT COUNTS   = " << distinct_counts << endl
	   << "MAX COUNT         = " << max_observed_count << endl
	   << "COUNTS OF 1       = " << counts_hist[1] << endl
	   << "MAX TERMS         = " << orig_max_terms << endl;
    
    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << endl;
    }
    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // BOOTSTRAPS

    if(bootstraps < 10)
      throw SMITHLABException("too few bootstraps, must be at least 10");

    if (VERBOSE) 
      cerr << "[BOOTSTRAP ESTIMATES]" << endl;
      
    vector<vector <double> > count_estimates;
    vector< vector<double> > sat_estimates;
    vector<double> lower_libsize, upper_libsize;
    estimates_bootstrap(VERBOSE, values,  bootstraps, orig_max_terms,
			diagonal, step_size, max_extrapolation, 
			max_val, val_step, count, count_estimates);
      
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // BOOTSTRAPS
    if (VERBOSE)
      cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
      
    vector<double> median_count_estimates;
    vector<double> count_upper_ci, count_lower_ci;
    return_median_and_ci(count_estimates, 1.0 - c_level, 
			 initial_count, median_count_estimates, 
			 count_lower_ci, count_upper_ci);

    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    if (VERBOSE) 
      cerr << "[WRITING OUTPUT]" << endl;
    
    write_predicted_curve(outfile, values_sum, c_level,
			  val_step, median_count_estimates,
			  count_lower_ci, count_upper_ci);
      
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
