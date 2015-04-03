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


int
main(const int argc, const char **argv) {

  try {

 
    bool VERBOSE = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;

    string outfile;

#ifdef HAVE_SAMTOOLS
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
#endif

    size_t max_num_points = 5;
    double tolerance = 1e-20;
    size_t max_iter = 0;
    size_t bootstraps = 100;
    double c_level = 0.95;
    double abundance_limit = 1e-8;

    bool UPPER_BOUND = false;

    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[1]),
                           "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("UPPER_BOUND", 'U', "compute upper bound, default "
		      "is lower bound", false, UPPER_BOUND);
    opt_parse.add_opt("abund_lim", 'a', "lower limit on abundance when computing "
		      "upper bounds (default: " + toa(abundance_limit) + ")",
		      false, abundance_limit);
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

    vector<string> leftover_args;
    opt_parse.parse(argc-1, argv+1, leftover_args);
    if (argc == 2 || opt_parse.help_requested()) {
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

    if(measure_moments.size() < 2*max_num_points)
      max_num_points = static_cast<size_t>(floor(measure_moments.size()/2));
    else
      measure_moments.resize(2*max_num_points);

    size_t n_points = ensure_pos_def_mom_seq(measure_moments, tolerance);

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

    double estimated_integral = 0.0;
    for(size_t i = 0; i < points.size(); i++)
      estimated_integral += counts_hist[1]*weights[i]/points[i];

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
