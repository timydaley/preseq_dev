/*    test_harris:
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

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::fixed;
using std::setprecision;


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
    size_t n_mixtures = 1;
    string alphas_string;
    string mus_string;
    string mixing_string;

    
    /* FLAGS */
    bool VERBOSE = false;
    //	bool SMOOTH_HISTOGRAM = false;	
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "",
			   "<sorted-bed-file>");
    opt_parse.add_opt("outfile", 'o', "output file for estimates",
		      false , outfile);
    opt_parse.add_opt("n_points",'p', "number of points for approximation",
		      false, num_points);
    opt_parse.add_opt("hist_max_terms",'h',"max terms in histogram",
		      false, hist_max_terms);
    opt_parse.add_opt("lib_size",'l', "library size",
		      false, lib_size);
    opt_parse.add_opt("n_mixtures", 'n', "number of mixtures",
		      true, n_mixtures);
    opt_parse.add_opt("mean", 'm', "mus, seperated by ,", true, mus_string);
    opt_parse.add_opt("alpha",'a',"alphas for NegBin dist",
    		      true, alphas_string);
    opt_parse.add_opt("mixing", 'x', "mixing parameters, must sum to 1",
		      true, mixing_string);
    opt_parse.add_opt("tol",'t',"numerical tolerance",
		      false, tolerance);
    opt_parse.add_opt("max_iter",'i',"maximum # iterations",
		      false, max_iter);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps to perform",
		      false, bootstraps);
    opt_parse.add_opt("CI",'c', "Confidence level",
		      false, CI);
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

    /*
    cerr << "READING IN PARAMS" << endl;
    vector<double> mixing;
    vector<double> distro_alphas;
    vector<double> mus;

    vector<string> alphas_parts = smithlab::split(alphas_string, ",");    
    vector<string> mus_parts  = smithlab::split(mus_string, ",");
    vector<string> mixing_parts = smithlab::split(mixing_string, ",");
    for(size_t i = 0; i < n_mixtures; i++){
      distro_alphas.push_back(atof(alphas_parts[i].c_str()));
      mus.push_back(atof(mus_parts[i].c_str()));
      mixing.push_back(atof(mixing_parts[i].c_str()));
    }    
    double mix_sum = accumulate(mixing.begin(), mixing.end(), 0.0);
    if( mix_sum!= 1){
      cerr << "Mixing parameters must sum to 1! \n";
      cerr << "Mixing sums to " << 
	mix_sum << ", need to reweight";
    }
    for(size_t i = 0; i < mixing.size(); i++)
      mixing[i] = mixing[i]/mix_sum;
    
    vector< vector<double> > alphas, betas;

    size_t iter = 0;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    for(size_t i = 0; i < num_points; i++)
      out << "alpha" << i << '\t';
    for(size_t i = 1; i < num_points; i++)
      out << "beta" << i << '\t';
    out << endl;


    while(iter < bootstraps){
      if(VERBOSE)
	cerr << "ITER " << iter << endl;
      iter++;
    // BUILD THE HISTOGRAM
    //    double mu = sampled_reads/lib_size;
      if(VERBOSE)
	cerr << "GENERATE SAMPLE" << endl; 
      vector<size_t> sample_counts(lib_size, 0);
      generate_NBD(mus, distro_alphas, mixing, sample_counts);
      const size_t max_observed_count = *std::max_element(sample_counts.begin(), sample_counts.end());
      vector<double> counts_hist(max_observed_count + 1, 0.0);
      for(size_t i = 0; i < sample_counts.size(); i++)
	counts_hist[sample_counts[i]]++;

      counts_hist[0] = 0;

      const double distinct_reads = accumulate(counts_hist.begin(), counts_hist.end(), 0.0);

      if (VERBOSE) {
	cerr << "LIBRARY_SIZE = " << lib_size << endl;
	cerr << "MU = ";
	for(size_t i = 0; i < mus.size(); i ++)
	  cerr << mus[i] << ", ";
	cerr << endl; 
	cerr << "ALPHA = ";
	for(size_t i = 0; i < distro_alphas.size(); i++)
	  cerr << distro_alphas[i] << ", ";
	cerr << endl; 
	cerr << "MIXING = ";
	for(size_t i = 0; i < mixing.size(); i++)
	  cerr << mixing[i] << ", ";
	cerr << endl;

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
	if(!isfinite(measure_moments.back())){
	  measure_moments.pop_back();
	  break;
	}
	indx++;
      }
      

      size_t n_points = std::min(num_points, static_cast<size_t>(floor(measure_moments.size()/2)));
      
      if(VERBOSE){
	cerr << "MOMENTS" << endl;
	for(size_t i = 0; i < measure_moments.size(); i++)
	  cerr << std::setprecision(16) << measure_moments[i] << endl;
	cerr << "OBSERVED_DISTINCT = " << accumulate(counts_hist.begin(), counts_hist.end(), 0.0) << endl;
      }
    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////

      MomentSequence quad_mom_seq(measure_moments);
      ZTNBD distro(1.0, 1.0);
      distro.EM_estim_params(tolerance, max_iter, counts_hist);
      vector<double> modified_moments;
      laguerre_modified_moments(measure_moments, distro.get_mu(),
				distro.get_alpha(),  num_points,
				modified_moments);
      const double k = 1/distro.get_alpha();
      const double p = distro.get_alpha()*distro.get_mu()/(1.0 + distro.get_alpha()*distro.get_mu());
      vector<double> modified_alpha, modified_beta;
      for(size_t i = 0; i < 2*num_points; i++)
	modified_alpha.push_back(p*(2*i + 1 + k));
      for(size_t i = 1; i < 2*num_points; i++)
	modified_beta.push_back(p*p*(i + k)*i);
      vector<double> full_alpha, full_beta;
      quad_mom_seq.modified_Chebyshev(VERBOSE, num_points, modified_alpha,
					   modified_beta, modified_moments,
					   full_alpha, full_beta);

      full_alpha.resize(num_points);
      full_beta.resize(num_points - 1);

      for(size_t i = 0; i < full_alpha.size(); i++)
	out << full_alpha[i] << '\t';
      for(size_t i = 0; i < full_beta.size(); i++)
	out << full_beta[i] << '\t';
      out << endl;
    }

    */


    
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
