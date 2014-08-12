/*    c_curve: plot a complexity curve by subsamping sequenced reads
 *    and counting UMIs
 *
 *    Copyright (C) 2012 University of Southern California and
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

#include <iomanip>
#include <numeric>
#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

#include "load_data_for_complexity.hpp"


using std::string;
using std::vector;
using std::endl;
using std::cerr;

using std::setw;
using std::fixed;
using std::setprecision;
using std::tr1::unordered_map;


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


int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    size_t lower_limit = 1000000;
    size_t upper_limit = 0;
    size_t step_size = 1000000;
    size_t mincount = 2;

    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;


    
#ifdef HAVE_SAMTOOLS
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
#endif


    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("c_curve", "plot a complexity curve by subsamping "
			   "sequenced reads and counting UMIs",
			   "<bed-file|bam-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("mincount",'m', "minimum # reads required to count a position",
		      false, mincount);
    opt_parse.add_opt("lower", 'l', "lower limit for samples", 
    		      false , lower_limit);
    opt_parse.add_opt("upper", 'u', "upper limit for samples", 
		      false , upper_limit);
    opt_parse.add_opt("step", 's', "step size for samples", 
    		      false , step_size);
    opt_parse.add_opt("verbose", 'v', "print more run information", 
		      false , VERBOSE);
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
    /**********************************************************************/
    
    // Setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    srand(time(0) + getpid());
    gsl_rng_set(rng, rand());
    
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


    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < counts_hist.size() &&
           counts_hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;

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
           << "COUNTS OF 1      = " << counts_hist[1] << endl;

    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
        if (counts_hist[i] > 0)
          cerr << i << '\t' << static_cast<size_t>(counts_hist[i]) << endl;
      cerr << endl;
    }

    //construct umi vector to sample from
    vector<size_t> umis;
    size_t umi = 1;
    for(size_t i = 1; i < counts_hist.size(); i++){
      for(size_t j = 0; j < counts_hist[i]; j++){
	for(size_t k = 0; k < i; k++)
	  umis.push_back(umi);
	umi++;
      }
    }
    assert(umis.size() == static_cast<size_t>(values_sum));
    
    if (upper_limit == 0)
      upper_limit = umis.size();
    else
      upper_limit = std::min(umis.size(), upper_limit);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    out << "total_reads" << "\t" << "distinct_reads" << endl;
    out << 0 << "\t" << 0 << endl;
    for (size_t i = lower_limit; i <= upper_limit; i += step_size) {
      if (VERBOSE)
	cerr << "sample size: " << i << endl;
      out << i << "\t" << sample_count_reads_w_mincount(rng, umis, 
							i, mincount) << endl;
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
