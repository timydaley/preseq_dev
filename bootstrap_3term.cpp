
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

void
resample_hist(const gsl_rng *rng,
              const vector<double> &counts_hist,
              vector<double> &sample_hist) {
    

  unsigned int distinct = 
    static_cast<unsigned int>(accumulate(counts_hist.begin(),
					 counts_hist.end(), 0.0));

  vector<unsigned int> sample_counts_hist(counts_hist.size(), 0);

        
  gsl_ran_multinomial(rng, counts_hist.size(), distinct,
		      &distinct_counts_hist.front(),
		      &sample_counts_hist.front());
    
  sample_hist.resize(counts_hist.size(), 0.0);
  for(size_t i = 0; i < sample_counts_hist.size(); i++)
    sample_hist[i] = static_cast<double>(sample_counts_hist[i]);
}


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
  const double phi = (1.0 + alpha*mu)/(alpha*mu);

  for(int n = 0; n < modified_moments.size(); n++){
    for(int l = 0; l <= n; l++){
      const double add_to_moment = 
	exp(gsl_sf_lngamma(n + k + 1) - gsl_sf_lngamma(n - l + 1)
	    - gsl_sf_lngamma(k + l + 1) + gsl_sf_lnfact(n)
	    - gsl_sf_lnfact(l) - (n - l)*log(phi)
	    + log(orig_moments[l]))*pow(-1, n + l);
      modified_moments[n] += add_to_moment; 
    } 
  } 
}

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
    double distro_alpha = 1.0;
    double distro_mu = 1.0;

    
    /* FLAGS */
    bool VERBOSE = false;
    //	bool SMOOTH_HISTOGRAM = false;	
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "", "<sorted-bed-file>");
    opt_parse.add_opt("outfile", 't', "output file name",
		      false , three_term_outfile);
    opt_parse.add_opt("n_points",'p', "number of points for approximation",
		      false, num_points);
    opt_parse.add_opt("lib_size",'l', "true library size",
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

      vector< vector<double> > alphas, betas;
      while(alphas.size() < bootstraps){
	vector<double> sample_hist;
	resample_hist(rng, counts_hist, sample_hist);

	vector<double> sample_moments;
  // mu_r = (r + 1)! n_{r+1} / n_1
	size_t indx = 1;
	while(sample_hist[indx] > 0  && indx <= sample_hist.size()){
	  sample_moments.push_back(exp(gsl_sf_lnfact(indx)
				       + log(sample_hist[indx])
				       - log(sample_hist[1])));
	  if(!std::isfinite(sample_moments.back())){
	    sample_moments.pop_back();
	    break;
	  }
	  indx++;
	}
      
	if(sample_moments.size()/2 >= n_points){

	  MomentSequence unmodCheb_mom_seq(sample_moments);
	  unmodCheb_mom_seq.unmodified_Chebyshev(VERBOSE);

	  if(unmodCheb_mom_seq.alpha.size() >= n_points){
	    unmod_mom_seq.alpha.resize(n_points);
	    unmod_mom_seq.beta.resize(n_points - 1);
	    alphas.push_back(unmod_mom_seq.alpha);
	    betas.push_back(unmod_mom_seq.beta);
	  }
	}
      }

      vector<double> median_alpha, median_beta;
      for(size_t j = 0; j < n_points; j++){
	vector<double> alpha_j;
	for(size_t i = 0; i < alphas.size(); i++)
	  alpha_j.push_back(alphas[i][j]);

	sort(alpha_j.begin(), alpha_j.end());
	median_alpha.push_back(gsl_stats_median_from_sorted_data(&alpha_j[0], 1, alpha_j.size()));
      }
      for(size_t j = 0; j < n_points - 1; j++){
	vector<double> beta_j;
	for(size_t i = 0; i < betas.size(); i++)
	  beta_j.push_back(betas[i][j]);

	sort(beta_j.begin(), beta_j.end());
	median_beta.push_back(gsl_stats_median_from_sorted_data(&beta_j[0], 1, beta_j.size()));
      }

      MomentSequence median_unmodCheb_mom_seq(median_alpha, median_beta);

      vector<double> points, weights;
      bool QUAD_SUCCESS = 
	median_unmodCheb_mom_seq.QR_quadrature_rules(VERBOSE, n_points,
						   tolerance, max_iter,
						   points, weights);

      double estimated_lib_size = 0.0;
      if(QUAD_SUCCESS){
	for(size_t i = 0; i < points.size(); i++)
	  estimated_lib_size += counts_hist[1]*weights[i]/points[i];
      }

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
      for(size_t i = 0 ; i < median_alpha.size(); i++)
	out << "median_alpha" << i << '\t';
      for(size_t i = 0; i < median_beta.size(); i++)
	out << "median_beta" << i << '\t';
      out << "estimated_lib_size" << endl;

      for(size_t i = 0; i < median_alpha.size(); i++)
	out << median_alpha[i] << '\t';
      for(size_t i = 0; i < median_beta.size() ; i++)
	out << median_beta[i] << '\t';
      out << estimated_lib_size << endl;
    
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
