/*    preseq: to predict properties of genomic sequencing libraries
 *
 *    Copyright (C) 2013 University of Southern California and
 *			 Andrew D. Smith and Timothy Daley
 *
 *    Authors: Timothy Daley, Victoria Helus, and Andrew Smith
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

#define PRESEQ_VERSION "0.1.0"

// AS: might not be good to depend on mapped read here
// TD: if we're including gc_extrap, we need the dependence
#include <MappedRead.hpp>

#include "continued_fraction.hpp"

using std::string;
using std::min;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::ifstream;

using std::setw;
using std::fixed;
using std::setprecision;
using std::priority_queue;
using std::tr1::unordered_map;


//////////////////////////////////////////////////////////////////////
// Data imputation
/////////////////////////////////////////////////////////////////////

/*
 * This code is used to deal with read data in BAM format.
 */
#ifdef HAVE_SAMTOOLS
// switching dependency on bamtools to samtools
#include <SAM.hpp>


const string mapper = "general";


static void
update_se_duplicate_counts_hist(const MappedRead &curr_mr,
                                const MappedRead &prev_mr,
                                const string input_file_name,
                                vector<double> &counts_hist,
                                size_t &current_count){
    // check if reads are sorted
    if (curr_mr.r.same_chrom(prev_mr.r) &&
        curr_mr.r.get_start() < prev_mr.r.get_start())
        throw SMITHLABException("locations unsorted in: " + input_file_name);
    
    if (!curr_mr.r.same_chrom(prev_mr.r) ||
        curr_mr.r.get_start() != prev_mr.r.get_start())
        // next read is new, update counts_hist to include current_count
    {
        // histogram is too small, resize
        if(counts_hist.size() < current_count + 1)
            counts_hist.resize(current_count + 1, 0.0);
        ++counts_hist[current_count];
        current_count = 1;
    }
    else // next read is same, update current_count
        ++current_count;
}


static size_t
load_counts_BAM_se(const string &input_file_name, vector<double> &counts_hist) {
    
    SAMReader sam_reader(input_file_name, mapper);
    if(!(sam_reader.is_good()))
        throw SMITHLABException("problem opening input file " + input_file_name);
    
    SAMRecord samr;
    sam_reader >> samr;
    size_t n_reads = 1;
    // resize vals_hist, make sure it starts out empty
    counts_hist.clear();
    counts_hist.resize(2, 0.0);
    size_t current_count = 1;
    
    MappedRead prev_mr, curr_mr;
    prev_mr = samr.mr;
    
    while (sam_reader >> samr, sam_reader.is_good())
    {
        if(samr.is_primary && samr.is_mapped)
            // only convert mapped and primary reads
        {
            // ignore unmapped reads & secondary alignments
            if(!(samr.is_mapping_paired) ||
               (samr.is_mapping_paired && samr.is_Trich)){
                //only count unpaired reads or the left mate of paired reads
                
                curr_mr = samr.mr;
                update_se_duplicate_counts_hist(curr_mr, prev_mr, input_file_name,
                                                counts_hist, current_count);
                
                // update number of reads and prev read
                ++n_reads;
                prev_mr = samr.mr;
            }
        }
    }
    
    // to account for the last read compared to the one before it.
    if(counts_hist.size() < current_count + 1)
        counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    
    return n_reads;
}

/********Below are functions for merging pair-end reads********/

static void
merge_mates(const size_t suffix_len, const size_t range,
            const GenomicRegion &one, const GenomicRegion &two,
            GenomicRegion &merged, int &len) {
    
    assert(one.same_chrom(two));
    const size_t read_start = min(one.get_start(), two.get_start());
    const size_t read_end = max(one.get_end(), two.get_end());
    
    len = read_end - read_start;
    
    if(len < 0){
        cerr << one << endl;
        cerr << two << endl;
        throw SMITHLABException("error merging reads");
    }
    
    merged = one;
    merged.set_start(read_start);
    merged.set_end(read_end);
    merged.set_score(one.get_score() + two.get_score());
    
    const string name(one.get_name());
    merged.set_name("FRAG:" + name.substr(0, name.size() - suffix_len));
    
}

// check if reads have same name & chrom
inline static bool
same_read(const size_t suffix_len, 
	  const MappedRead &a, const MappedRead &b) {
  const string sa(a.r.get_name());
  const string sb(b.r.get_name());
  return (std::equal(sa.begin(), sa.end() - suffix_len, sb.begin())
	  && a.r.same_chrom(b.r));
}


static bool
update_pe_duplicate_counts_hist(const GenomicRegion &curr_gr,
                                const GenomicRegion &prev_gr,
                                vector<double> &counts_hist,
                                size_t &current_count){
    // check if reads are sorted
    if (curr_gr.same_chrom(prev_gr) &&
        curr_gr.get_start() < prev_gr.get_start()
        && curr_gr.get_end() < prev_gr.get_end()){
        return false;
    }
    if (!curr_gr.same_chrom(prev_gr) ||
        curr_gr.get_start() != prev_gr.get_start() ||
        curr_gr.get_end() != prev_gr.get_end())
        // next read is new, update counts_hist to include current_count
    {
        // histogram is too small, resize
        if(counts_hist.size() < current_count + 1)
            counts_hist.resize(current_count + 1, 0.0);
        ++counts_hist[current_count];
        current_count = 1;
    }
    else // next read is same, update current_count
        ++current_count;
    
    return true;
}

/////comparison function for priority queue/////////////////

/**************** FOR CLARITY BELOW WHEN COMPARING READS *************/
static inline bool
chrom_greater(const GenomicRegion &a, const GenomicRegion &b) {
    return a.get_chrom() > b.get_chrom();
}
static inline bool
same_start(const GenomicRegion &a, const GenomicRegion &b) {
    return a.get_start() == b.get_start();
}
static inline bool
start_greater(const GenomicRegion &a, const GenomicRegion &b) {
    return a.get_start() > b.get_start();
}
static inline bool
end_greater(const GenomicRegion &a, const GenomicRegion &b) {
    return a.get_end() > b.get_end();
}
/******************************************************************************/


struct GenomicRegionOrderChecker {
    bool operator()(const GenomicRegion &prev, const GenomicRegion &gr) const {
        return start_check(prev, gr);
    }
    static bool
    is_ready(const priority_queue<GenomicRegion, vector<GenomicRegion>, GenomicRegionOrderChecker> &pq,
             const GenomicRegion &gr, const size_t max_width) {
        return !pq.top().same_chrom(gr) || pq.top().get_end() + max_width < gr.get_start();
    }
    static bool
    start_check(const GenomicRegion &prev, const GenomicRegion &gr) {
        return (chrom_greater(prev, gr)
                || (prev.same_chrom(gr) && start_greater(prev, gr))
                || (prev.same_chrom(gr) && same_start(prev, gr) && end_greater(prev, gr)));
    }
};



static void empty_pq(GenomicRegion &curr_gr, GenomicRegion &prev_gr,
                     size_t &current_count, vector<double> &counts_hist,
                     priority_queue<GenomicRegion, vector<GenomicRegion>,
                                    GenomicRegionOrderChecker> &read_pq,
                     const string &input_file_name ){
    
    curr_gr = read_pq.top();
    //	       cerr << "outputting from queue : " << read_pq.top() << endl;
    read_pq.pop();
    
    //update counts hist
    bool UPDATE_SUCCESS
    = update_pe_duplicate_counts_hist(curr_gr, prev_gr, counts_hist,
                                      current_count);
    if(!UPDATE_SUCCESS){
        //cerr << "prev = " << prev_gr << endl;
        //cerr << "curr = " << curr_gr << endl;
        //cerr << "priority queue : " << endl;
        //while(	 !(read_pq.empty()) ){
          //  cerr << read_pq.top() << endl;
          //  read_pq.pop();
        //}
        throw SMITHLABException("reads unsorted in " + input_file_name);
    }
    prev_gr = curr_gr;
}



static size_t
load_counts_BAM_pe(const bool VERBOSE,
                   const string &input_file_name,
                   const size_t MAX_SEGMENT_LENGTH,
                   const size_t MAX_READS_TO_HOLD,
                   vector<double> &counts_hist) {
    
    SAMReader sam_reader(input_file_name, mapper);
    // check sam_reader
    if(!(sam_reader.is_good()))
        throw SMITHLABException("problem opening input file " + input_file_name);
    
    SAMRecord samr;
    size_t n_reads = 0;
    // resize vals_hist, make sure it starts out empty
    counts_hist.clear();
    counts_hist.resize(2, 0.0);
    size_t current_count = 1;
    size_t suffix_len = 0;
    size_t n_merged = 0;
    size_t n_unpaired = 0;
    
    GenomicRegion curr_gr, prev_gr;
    
    std::priority_queue<GenomicRegion, vector<GenomicRegion>,
    GenomicRegionOrderChecker> read_pq;
    
    unordered_map<string, SAMRecord> dangling_mates;
    
    while ((sam_reader >> samr, sam_reader.is_good()))
    {
      if(samr.is_primary && samr.is_mapped){
	// only convert mapped and primary reads
	if (samr.is_mapping_paired){
	  const string read_name
	    = samr.mr.r.get_name().substr(0, samr.mr.r.get_name().size() - suffix_len);
                
	  if (dangling_mates.find(read_name) != dangling_mates.end()){
	    // other end is in dangling mates, merge the two mates
	    if(same_read(suffix_len, samr.mr, dangling_mates[read_name].mr)){
	      if (samr.is_Trich) std::swap(samr, dangling_mates[read_name]);
                    
	      GenomicRegion merged;
	      int len = 0;
	      merge_mates(suffix_len, MAX_SEGMENT_LENGTH,
			  dangling_mates[read_name].mr.r, samr.mr.r, merged, len);
	    // merge success!
	      if (len >= 0 && len <= static_cast<int>(MAX_SEGMENT_LENGTH)){
	      // first iteration
		if(n_reads == 0){
		  prev_gr = merged;
		  ++n_reads;
		  ++n_merged;
		}
		else{
		  ++n_reads;
		  ++n_merged;
		  read_pq.push(merged);
                            
		  if(!(read_pq.empty()) &&
		     GenomicRegionOrderChecker::is_ready(read_pq, merged, MAX_SEGMENT_LENGTH)) {
		  //begin emptying priority queue
		    while(!(read_pq.empty()) &&
			  GenomicRegionOrderChecker::is_ready(read_pq, merged,
							      MAX_SEGMENT_LENGTH) ){
		      empty_pq(curr_gr, prev_gr, current_count,
			       counts_hist, read_pq, input_file_name);
		    }//end while loop
		  }//end statement for emptying priority queue
		}
		dangling_mates.erase(read_name);
	     

	      }//end if statement for if merge is successful
	      else{
		cerr << "problem with read " << read_name << endl;

		throw SMITHLABException("merge unsuccessful");
	      }
	    }
	    else{
	      // problem mergin reads from "different" chrom like chr1 & chr1_gl000191_random
	      // flagged as proper pair, but not
	      read_pq.push(samr.mr.r);
	      read_pq.push(dangling_mates[read_name].mr.r);
	    }
	  }//end if statement for if read is in dangling mates
	  else{	// other end not in dangling mates, add this read to dangling mates.
	    dangling_mates[read_name] = samr;
	  }
	}
	else{ // read is unpaired, put in queue
	  if(n_reads == 0){
	    ++n_reads;
	    ++n_unpaired;
	    prev_gr = samr.mr.r;
	  }
	  else{ 
	    ++n_reads;
	    ++n_unpaired;
	    read_pq.push(samr.mr.r);
 
	    if(!(read_pq.empty()) &&
	       GenomicRegionOrderChecker::is_ready(read_pq, samr.mr.r, MAX_SEGMENT_LENGTH)) {
	     
	      while(!(read_pq.empty()) &&
		    GenomicRegionOrderChecker::is_ready(read_pq, samr.mr.r,
							MAX_SEGMENT_LENGTH) ){
                              
		empty_pq(curr_gr, prev_gr, current_count,
			 counts_hist, read_pq, input_file_name);
	      }
                    
	    }
                
	  }

	} 
            
	// dangling mates is too large, flush dangling_mates of reads
	// on different chroms and too far away
	if (dangling_mates.size() > MAX_READS_TO_HOLD){
	  if(VERBOSE)
	    cerr << "dangling mates too large, emptying" << endl;
                
	  unordered_map<string, SAMRecord> tmp;
	  for (unordered_map<string, SAMRecord>::iterator
		 itr = dangling_mates.begin();
	       itr != dangling_mates.end(); ++itr)
	    if (itr->second.mr.r.get_chrom() != samr.mr.r.get_chrom()
		|| (itr->second.mr.r.get_chrom() == samr.mr.r.get_chrom()
		    && itr->second.mr.r.get_end() + MAX_SEGMENT_LENGTH <
		    samr.mr.r.get_start())) 
	      {
		if(itr->second.seg_len >= 0)
		  read_pq.push(itr->second.mr.r);
	      }
	    else
	      tmp[itr->first] = itr->second;
                
	  std::swap(tmp, dangling_mates);
	}
            
      }
    }
    
    // empty dangling mates of any excess reads
    while (!dangling_mates.empty()) {
      read_pq.push(dangling_mates.begin()->second.mr.r);
      dangling_mates.erase(dangling_mates.begin());
    }
  
    //final iteration
    while(!read_pq.empty()){
      empty_pq(curr_gr, prev_gr, current_count,
	       counts_hist, read_pq, input_file_name);
    }
    
    if(counts_hist.size() < current_count + 1)
      counts_hist.resize(current_count + 1, 0.0);

    ++counts_hist[current_count];

    assert((read_pq.empty()));
    
    return n_reads;
}

#endif


/*
 this code is for BED file input
 */

static void
update_se_duplicate_counts_hist(const GenomicRegion &curr_gr,
                                const GenomicRegion &prev_gr,
                                const string input_file_name,
                                vector<double> &counts_hist,
                                size_t &current_count){
    // check if reads are sorted
    if (curr_gr.same_chrom(prev_gr) &&
        curr_gr.get_start() < prev_gr.get_start())
        throw SMITHLABException("locations unsorted in: " + input_file_name);
    
    if (!curr_gr.same_chrom(prev_gr) ||
        curr_gr.get_start() != prev_gr.get_start())
        // next read is new, update counts_hist to include current_count
    {
        // histogram is too small, resize
        if(counts_hist.size() < current_count + 1)
            counts_hist.resize(current_count + 1, 0.0);
        ++counts_hist[current_count];
        current_count = 1;
    }
    else // next read is same, update current_count
        ++current_count;
}

static size_t
load_counts_BED_se(const string input_file_name, vector<double> &counts_hist) {
    // resize vals_hist
    counts_hist.clear();
    counts_hist.resize(2, 0.0);
    
    std::ifstream in(input_file_name.c_str());
    if (!in)
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    GenomicRegion curr_gr, prev_gr;
    if (!(in >> prev_gr))
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    size_t n_reads = 1;
    size_t current_count = 1;
    
    while (in >> curr_gr) {
        update_se_duplicate_counts_hist(curr_gr, prev_gr, input_file_name,
                                        counts_hist, current_count);
        ++n_reads;
        prev_gr.swap(curr_gr);
    }
    
    
    // to account for the last read compared to the one before it.
    if(counts_hist.size() < current_count + 1)
        counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    
    return n_reads;
}


static size_t
load_counts_BED_pe(const string input_file_name, vector<double> &counts_hist) {
    
    // resize vals_hist
    counts_hist.clear();
    counts_hist.resize(2, 0.0);
    
    std::ifstream in(input_file_name.c_str());
    if (!in)
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    GenomicRegion curr_gr, prev_gr;
    if (!(in >> prev_gr))
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    size_t n_reads = 1;
    size_t current_count = 1;
    
    //read in file and compare each gr with the one before it
    while (in >> curr_gr) {
        bool UPDATE_SUCCESS =
        update_pe_duplicate_counts_hist(curr_gr, prev_gr,
                                        counts_hist, current_count);
        if(!UPDATE_SUCCESS){
           // cerr << "prev = " << prev_gr << endl;
           // cerr << "curr = " << curr_gr << endl;
            throw SMITHLABException("reads unsorted in " + input_file_name);
        }
        
        ++n_reads;
        prev_gr.swap(curr_gr);
    }
    
    if(counts_hist.size() < current_count + 1)
        counts_hist.resize(current_count + 1, 0.0);
    
    // to account for the last read compared to the one before it.
    ++counts_hist[current_count];
    
    
    return n_reads;
    
    
}

/*
 text file input
 */

static size_t
load_counts(const string input_file_name, vector<double> &counts_hist) {
    
    std::ifstream in(input_file_name.c_str());
    if (!in) // if file doesn't open
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    size_t n_reads = 0;
    while(!in.eof()){
        string buffer;
        getline(in, buffer);
        
        std::istringstream iss(buffer);
        if(iss.good()){
            double val;
            iss >> val;
            if(val > 0) {
                const size_t count = static_cast<size_t>(val);
                // histogram is too small, resize
                if(counts_hist.size() < count + 1)
                    counts_hist.resize(count + 1, 0.0);
                ++counts_hist[count];
                n_reads += count;
            }
            else if(val != 0)
                throw SMITHLABException("problem reading file at line " + (n_reads + 1));
        }
        in.peek();
    }
    in.close();
    
    return n_reads;
}

//returns number of reads from file containing counts histogram
static size_t
load_histogram(const string &filename, vector<double> &counts_hist) {
    
    counts_hist.clear();
    
    std::ifstream in(filename.c_str());
    if (!in) //if file doesn't open
        throw SMITHLABException("could not open histogram: " + filename);
    
    size_t n_reads = 0;
    size_t line_count = 0ul, prev_read_count = 0ul;
    string buffer;
    while (getline(in, buffer)) {
        ++line_count;
        size_t read_count = 0ul;
        double frequency = 0.0;
        std::istringstream is(buffer);
        // error reading input
        if (!(is >> read_count >> frequency))
            throw SMITHLABException("bad histogram line format:\n" +
                                    buffer + "\n(line " + toa(line_count) + ")");
        // histogram is out of order
        if (read_count < prev_read_count)
            throw SMITHLABException("bad line order in file " +
                                    filename + "\n(line " +
                                    toa(line_count) + ")");
        counts_hist.resize(read_count + 1, 0.0);
        counts_hist[read_count] = frequency;
        prev_read_count = read_count;
        n_reads += static_cast<size_t>(read_count*frequency);
    }
    
    return n_reads;
}

/////////////////////////////////////////////////////////
// Confidence interval stuff

static inline double
alpha_log_confint_multiplier(const double estimate,
                             const double variance, const double alpha) {
    const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
    return exp(inv_norm_alpha*
               sqrt(log(1.0 + variance/pow(estimate, 2))));
}


/*static void
ci_given_estimates(const vector<vector<double> > &bootstrap_estimates,
                   const double ci_level, const vector<double> &yield_estimates,
                   vector<double> &lower_ci_lognormal,
                   vector<double> &upper_ci_lognormal) {
    
    lower_ci_lognormal.clear();
    upper_ci_lognormal.clear();
    
    const double alpha = 1.0 - ci_level;
    assert(!bootstrap_estimates.empty());
    
    const size_t n_est = bootstrap_estimates.size();
    vector<double> estimates_row(bootstrap_estimates.size(), 0.0);
    for (size_t i = 0; i < bootstrap_estimates[0].size(); i++) {
        
        // estimates is in wrong order, work locally
        for (size_t k = 0; k < n_est; ++k)
            estimates_row[k] = bootstrap_estimates[k][i];
        
        // sort to get confidence interval
        const double variance = gsl_stats_variance(&estimates_row[0], 1, n_est);
        const double confint_mltr =
        alpha_log_confint_multiplier(yield_estimates[i], variance, alpha);
        lower_ci_lognormal.push_back(yield_estimates[i]/confint_mltr);
        upper_ci_lognormal.push_back(yield_estimates[i]*confint_mltr);
    }
}*/

static void
median_and_ci(const vector<double> &estimates,
              const double ci_level,
              double &median_estimate,
              double &lower_ci_estimate,
              double &upper_ci_estimate){
    assert(!estimates.empty());
    const double alpha = 1.0 - ci_level;
    const size_t n_est = estimates.size();
    vector<double> sorted_estimates(estimates);
    sort(sorted_estimates.begin(), sorted_estimates.end());
    median_estimate =
    gsl_stats_median_from_sorted_data(&sorted_estimates[0], 1, n_est);
    const double variance = gsl_stats_variance(&sorted_estimates[0], 1, n_est);
    const double confint_mltr =
    alpha_log_confint_multiplier(median_estimate, variance, alpha);
    
    lower_ci_estimate = median_estimate/confint_mltr;
    upper_ci_estimate = median_estimate*confint_mltr;
    
}

static void
vector_median_and_ci(const vector<vector<double> > &bootstrap_estimates,
                     const double ci_level, vector<double> &yield_estimates,
                     vector<double> &lower_ci_lognormal,
                     vector<double> &upper_ci_lognormal) {
    
    yield_estimates.clear();
    lower_ci_lognormal.clear();
    upper_ci_lognormal.clear();
    assert(!bootstrap_estimates.empty());
    
    const size_t n_est = bootstrap_estimates.size();
    vector<double> estimates_row(bootstrap_estimates.size(), 0.0);
    for (size_t i = 0; i < bootstrap_estimates[0].size(); i++) {
        
        // estimates is in wrong order, work locally on const val
        for (size_t k = 0; k < n_est; ++k)
            estimates_row[k] = bootstrap_estimates[k][i];
        
        double median_estimate, lower_ci_estimate, upper_ci_estimate;
        median_and_ci(estimates_row, ci_level, median_estimate,
                      lower_ci_estimate, upper_ci_estimate);
        sort(estimates_row.begin(), estimates_row.end());
        
        yield_estimates.push_back(median_estimate);
        lower_ci_lognormal.push_back(lower_ci_estimate);
        upper_ci_lognormal.push_back(upper_ci_estimate);
    }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
/////  LC_EXTRAP MODE BELOW HERE
/////



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
sample_count_distinct(const gsl_rng *rng,
                      const vector<size_t> &full_umis,
                      const size_t sample_size) {
    vector<size_t> sample_umis(sample_size);
    gsl_ran_choose(rng, (size_t *)&sample_umis.front(), sample_size,
                   (size_t *)&full_umis.front(), full_umis.size(),
                   sizeof(size_t));
    double count = 1.0;
    for (size_t i = 1; i < sample_umis.size(); i++)
        if(sample_umis[i] != sample_umis[i-1])
            count++;
    
    return count;
}


// check if estimates are finite, increasing, and concave
static bool
check_yield_estimates(const vector<double> &estimates) {
    
    if (estimates.empty())
        return false;
    
    // make sure that the estimate is increasing in the time_step and is
    // below the initial distinct per step_size
    if (!finite(accumulate(estimates.begin(), estimates.end(), 0.0)))
        return false;
    
    for (size_t i = 1; i < estimates.size(); ++i)
        if ((estimates[i] < estimates[i - 1]) ||
            (i >= 2 && (estimates[i] - estimates[i - 1] >
                        estimates[i - 1] - estimates[i - 2])) ||
            (estimates[i] < 0.0))
            return false;
    
    return true;
}


void
lc_extrap_bootstrap(const bool VERBOSE, const vector<double> &orig_hist,
                    const size_t bootstraps, const size_t orig_max_terms,
                    const int diagonal, const double step_size,
                    const double max_extrapolation, const size_t max_iter,
                    vector< vector<double> > &bootstrap_estimates) {
    // clear returning vectors
    bootstrap_estimates.clear();
    
    //setup rng
    srand(time(0) + getpid());
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, rand());
    
    double vals_sum = 0.0;
    for(size_t i = 0; i < orig_hist.size(); i++)
        vals_sum += orig_hist[i]*i;
    
    const double initial_distinct = accumulate(orig_hist.begin(), orig_hist.end(), 0.0);
    
    
    vector<size_t> orig_hist_distinct_counts;
    vector<double> distinct_orig_hist;
    for(size_t i = 0; i < orig_hist.size(); i++){
        if(orig_hist[i] > 0){
            orig_hist_distinct_counts.push_back(i);
            distinct_orig_hist.push_back(orig_hist[i]);
        }
    }
    
    for (size_t iter = 0;
         (iter < max_iter && bootstrap_estimates.size() < bootstraps);
         ++iter) {
        
        vector<double> yield_vector;
        vector<double> hist;
        resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);
        
        double sample_vals_sum = 0.0;
        for(size_t i = 0; i < hist.size(); i++)
            sample_vals_sum += i*hist[i];
        
        // const double sample_max_val = max_extrapolation/sample_vals_sum;
        //   const double sample_val_step = step_size/sample_vals_sum;
        
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
        assert(umis.size() == static_cast<size_t>(sample_vals_sum));
        
        // interpolate complexity curve by random sampling w/out replacement
        size_t upper_limit = static_cast<size_t>(sample_vals_sum);
        size_t step = static_cast<size_t>(step_size);
        size_t sample = step;
        while(sample < upper_limit){
            yield_vector.push_back(sample_count_distinct(rng, umis, sample));
            sample += step;
        }
        
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
        lower_cfa(diagonal, max_terms, step_size, max_extrapolation);
        
        const ContinuedFraction
        lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));
        
        //extrapolate the curve start
        if (lower_cf.is_valid()){
            double sample_size = static_cast<double>(sample);
            while(sample_size < max_extrapolation){
	      const double one_minus_fold_extrap = (sample_size - sample_vals_sum)/sample_vals_sum;
	      assert(one_minus_fold_extrap >= 0.0);
	      yield_vector.push_back(initial_distinct 
				     + one_minus_fold_extrap*lower_cf(one_minus_fold_extrap));
	      sample_size += step_size;
            }
            
            // SANITY CHECK
            if (check_yield_estimates(yield_vector)) {
                bootstrap_estimates.push_back(yield_vector);
                if (VERBOSE) cerr << lower_cf.degree << ", ";
            }
            else if (VERBOSE){
                cerr << "_";
            }
        }
        else if (VERBOSE){
            cerr << "_";
        }
        
    }
    if (VERBOSE)
        cerr << endl;
    if (bootstrap_estimates.size() < bootstraps)
        throw SMITHLABException("too many iterations, poor sample");
}

static bool
lc_extrap_single_estimate(const bool VERBOSE, vector<double> &hist,
                          size_t max_terms, const int diagonal,
                          const double step_size, const double max_extrapolation,
                          vector<double> &yield_estimate) {
    
    //setup rng
    srand(time(0) + getpid());
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, rand());
    
    
    yield_estimate.clear();
    double vals_sum = 0.0;
    for(size_t i = 0; i < hist.size(); i++)
        vals_sum += i*hist[i];
    const double initial_distinct = accumulate(hist.begin(), hist.end(), 0.0);
    
    // const double max_val = max_extrapolation/vals_sum;
    // const double val_step = step_size/vals_sum;
    
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
    while(sample < upper_limit){
        yield_estimate.push_back(sample_count_distinct(rng, umis, sample));
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
    lower_cfa(diagonal, max_terms, step_size, max_extrapolation);
    
    const ContinuedFraction
    lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));
    
    // extrapolate curve
    if (lower_cf.is_valid()){
        double sample_size = static_cast<double>(sample);
        while(sample_size < max_extrapolation){
            const double one_minus_fold_extrap = (sample_size - vals_sum)/vals_sum;
            assert(one_minus_fold_extrap >= 0.0);
            yield_estimate.push_back(initial_distinct + one_minus_fold_extrap*lower_cf(one_minus_fold_extrap));
            sample_size += step_size;
        }
    }
    else{
        // FAIL!
        // lower_cf unacceptable, need to bootstrap to obtain estimates
        return false;
    }
    
    if (VERBOSE) {
        cerr << "CF_OFFSET_COEFF_ESTIMATES" << endl;
        copy(lower_cf.offset_coeffs.begin(), lower_cf.offset_coeffs.end(),
             std::ostream_iterator<double>(cerr, "\n"));
        
        cerr << "CF_COEFF_ESTIMATES" << endl;
        copy(lower_cf.cf_coeffs.begin(), lower_cf.cf_coeffs.end(),
             std::ostream_iterator<double>(cerr, "\n"));
    }
    
    // SUCCESS!!
    return true;
}

static double
GoodToulmin2xExtrap(const vector<double> &counts_hist){
    double two_fold_extrap = 0.0;
    for(size_t i = 0; i < counts_hist.size(); i++)
        two_fold_extrap += pow(-1.0, i + 1)*counts_hist[i];
    
    return two_fold_extrap;
}


static void
write_predicted_complexity_curve(const string outfile,
                                 const double c_level, const double step_size,
                                 const vector<double> &yield_estimates,
                                 const vector<double> &yield_lower_ci_lognormal,
                                 const vector<double> &yield_upper_ci_lognormal) {
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    out << "TOTAL_READS\tEXPECTED_DISTINCT\t"
    << "LOGNORMAL_LOWER_" << 100*c_level << "%CI\t"
    << "LOGNORMAL_UPPER_" << 100*c_level << "%CI" << endl;
    
    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    out.precision(1);
    
    out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
    for (size_t i = 0; i < yield_estimates.size(); ++i)
        out << (i + 1)*step_size << '\t'
        << yield_estimates[i] << '\t'
        << yield_lower_ci_lognormal[i] << '\t' << yield_upper_ci_lognormal[i] << endl;
}




static void
lc_extrap(const bool VERBOSE,
          const bool VALS_INPUT,
          const bool PAIRED_END,
          const bool HIST_INPUT,
          const bool SINGLE_ESTIMATE,
#ifdef HAVE_SAMTOOLS
          const bool BAM_FORMAT_INPUT,
#endif
          const size_t MIN_REQUIRED_COUNTS,
          size_t orig_max_terms,
          const double max_extrapolation,
          double step_size,
          const size_t bootstraps,
          const int diagonal,
          const double c_level,
          const string &input_file_name,
          const string &outfile,
          const size_t MAX_SEGMENT_LENGTH) {
    
    vector<double> counts_hist;
    size_t n_reads = 0;
    
    // LOAD VALUES
    if(HIST_INPUT){
        if(VERBOSE)
            cerr << "INPUT_HIST" << endl;
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
        const size_t MAX_READS_TO_HOLD = 100000;
        n_reads = load_counts_BAM_pe(VERBOSE, input_file_name, MAX_SEGMENT_LENGTH,
                                     MAX_READS_TO_HOLD, counts_hist);
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
        step_size = std::max(step_size, step_size*round(n_reads/(20*step_size)));
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
        throw SMITHLABException("Library expected to saturate in doubling of size, unable to extrapolate");
    
    
    size_t total_reads = 0;
    for(size_t i = 0; i < counts_hist.size(); i++){
        total_reads += i*counts_hist[i];
      //  cerr << "total reads " << total_reads << endl;
    }
    //assert(total_reads == n_reads);
    
    
    
    // catch if all reads are distinct or sample sufficiently deep
    if (max_observed_count < MIN_REQUIRED_COUNTS)
        throw SMITHLABException("sample not sufficiently deep or duplicates removed");
    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // ESTIMATE COMPLEXITY CURVE
    
    if(VERBOSE)
        cerr << "[ESTIMATING YIELD CURVE]" << endl;
    vector<double> yield_estimates;
 
    
    if(SINGLE_ESTIMATE){
      bool SINGLE_ESTIMATE_SUCCESS =
	lc_extrap_single_estimate(VERBOSE, counts_hist, orig_max_terms, diagonal,
				  step_size, max_extrapolation, yield_estimates);
        // IF FAILURE, EXIT
        if(!SINGLE_ESTIMATE_SUCCESS)
            throw SMITHLABException("SINGLE ESTIMATE FAILED, NEED TO RUN FULL MODE FOR ESTIMATES");
        
        std::ofstream of;
        if (!outfile.empty()) of.open(outfile.c_str());
        std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
        
        out << "TOTAL_READS\tEXPECTED_DISTINCT" << endl;
        
        out.setf(std::ios_base::fixed, std::ios_base::floatfield);
        out.precision(1);
        
        out << 0 << '\t' << 0 << endl;
        for (size_t i = 0; i < yield_estimates.size(); ++i)
            out << (i + 1)*step_size << '\t'
            << yield_estimates[i] << endl;
        
    }
    else{
        if (VERBOSE)
            cerr << "[BOOTSTRAPPING HISTOGRAM]" << endl;
        
        const size_t max_iter = 4*bootstraps;
        
        vector<vector <double> > bootstrap_estimates;
        lc_extrap_bootstrap(VERBOSE, counts_hist, bootstraps, orig_max_terms,
                            diagonal, step_size, max_extrapolation, max_iter,
                            bootstrap_estimates);
        
        
        /////////////////////////////////////////////////////////////////////
        if (VERBOSE)
            cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
        
        // yield ci
        vector<double> yield_upper_ci_lognormal, yield_lower_ci_lognormal;
	     
        //if(!SINGLE_ESTIMATE_SUCCESS){
            // use bootstrap estimates to obtain median estimates
            vector_median_and_ci(bootstrap_estimates, c_level, yield_estimates,
                                 yield_lower_ci_lognormal, yield_upper_ci_lognormal);
	    /*    }
        else{
	
            // use single estimates as the expected complexity curve
            ci_given_estimates(bootstrap_estimates, c_level, yield_estimates,
                               yield_lower_ci_lognormal, yield_upper_ci_lognormal);
			       }*/
        
        /////////////////////////////////////////////////////////////////////
        if (VERBOSE)
            cerr << "[WRITING OUTPUT]" << endl;
        
        write_predicted_complexity_curve(outfile, c_level, step_size,
                                         yield_estimates, yield_lower_ci_lognormal,
                                         yield_upper_ci_lognormal);
    }
    
}




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
/////  C_CURVE BELOW HERE
/////



static void c_curve (const bool VERBOSE,
                     const bool VALS_INPUT,
                     const bool PAIRED_END,
                     const bool HIST_INPUT,
#ifdef HAVE_SAMTOOLS
                     const bool BAM_FORMAT_INPUT,
#endif
                     double step_size,
                     size_t upper_limit,
                     const string &input_file_name,
                     const string &outfile,
                     const size_t MAX_SEGMENT_LENGTH){
    
    
    // Setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default); // use default type
    srand(time(0) + getpid()); //give the random fxn a new seed
    gsl_rng_set(rng, rand()); //initialize random number generator with the seed
    
    vector<double> counts_hist;
    size_t n_reads = 0;
    
    // LOAD VALUES
    if(HIST_INPUT){
        if(VERBOSE)
            cerr << "INPUT_HIST" << endl;
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
        const size_t MAX_READS_TO_HOLD = 100000;
        n_reads = load_counts_BAM_pe(VERBOSE, input_file_name, MAX_SEGMENT_LENGTH,
                                     MAX_READS_TO_HOLD, counts_hist);
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
    
    
    size_t total_reads = 0;
    for(size_t i = 0; i < counts_hist.size(); i++)
        total_reads += i*counts_hist[i];
    
    const size_t distinct_counts =
    static_cast<size_t>(std::count_if(counts_hist.begin(), counts_hist.end(),
                                      bind2nd(std::greater<double>(), 0.0)));
    if (VERBOSE)
        cerr << "TOTAL READS     = " << n_reads << endl
        << "DISTINCT READS  = " << distinct_reads << endl
        << "DISTINCT COUNTS = " << distinct_counts << endl
        << "MAX COUNT       = " << max_observed_count << endl
        << "COUNTS OF 1     = " << counts_hist[1] << endl;
    
    if (VERBOSE) {
        // OUTPUT THE ORIGINAL HISTOGRAM
        cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
        for (size_t i = 0; i < counts_hist.size(); i++)
            if (counts_hist[i] > 0)
                cerr << i << '\t' << static_cast<size_t>(counts_hist[i]) << endl;
        cerr << endl;
    }
    
    //if(total_reads != n_reads){
    //    cerr << "total reads = " << total_reads << endl;
    //    cerr << "n_reads     = " << n_reads << endl;
    //}
    
    //assert(total_reads == n_reads);
    
    
    //construct umi vector to sample from
    vector<size_t> full_umis;
    size_t umi = 1;
    for(size_t i = 1; i < counts_hist.size(); i++){
        for(size_t j = 0; j < counts_hist[i]; j++){
            for(size_t k = 0; k < i; k++)
                full_umis.push_back(umi);
            umi++;
        }
    }
    
    assert(full_umis.size() == n_reads);
    
    if (upper_limit == 0)
        upper_limit = n_reads; //set upper limit to equal the number of molecules
    
    
	//handles output of c_curve
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    
	//prints the complexity curve
    out << "total_reads" << "\t" << "distinct_reads" << endl;
    out << 0 << '\t' << 0 << endl;
    for (size_t i = step_size; i <= upper_limit; i += step_size) {
        if (VERBOSE)
            cerr << "sample size: " << i << endl;
        out << i << "\t" << sample_count_distinct(rng, full_umis, i) << endl;
    }
}






static int usage()
{
    cerr << "\nProgram: preseq (applications for analyzing library complexity)" << endl;
    cerr << "Version: "<< PRESEQ_VERSION << "\n\n";
    cerr << "Usage: preseq <command> [OPTIONS]\n\n";
    cerr << "Command: c_curve          generate complexity curve for a library\n";
    cerr << "         lc_extrap        predict the yield for future experiments\n";
    // cerr << "         gc_extrap        extrapolate genomic complexity";
    cerr << "\n\n";
    return 0;
    //   string input_file_name2 = "SRR726645.sort.mr";
    
}
size_t upper_limit = 0;



int
main(const int argc, const char **argv) {
    
    try {
        
        const size_t MIN_REQUIRED_COUNTS = 4;
        
        /* FILES */
        string outfile;
        
        size_t orig_max_terms = 100;
        double max_extrapolation = 1.0e10;
        size_t upper_limit = 0;
        
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
        
        if (argc < 2) return usage();
        else if (strcmp(argv[1], "lc_extrap") == 0) {
            
            /********** GET COMMAND LINE ARGUMENTS  FOR LC EXTRAP ***********/
            OptionParser opt_parse(strip_path(argv[1]),
                                   "", "<sorted-bed-file>");
            opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
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
                              "quick mode, estimate yield without bootstrapping for confidence intervals",
                              false, SINGLE_ESTIMATE);
            
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
            /******************************************************************/
            
            if (argc > 2) {
                //    cerr << leftover_args.front();
                lc_extrap(VERBOSE,
                          VALS_INPUT,
                          PAIRED_END,
                          HIST_INPUT,
                          SINGLE_ESTIMATE,
#ifdef HAVE_SAMTOOLS
                          BAM_FORMAT_INPUT,
#endif
                          MIN_REQUIRED_COUNTS,
                          orig_max_terms,
                          max_extrapolation,
                          step_size,
                          bootstraps,
                          diagonal,
                          c_level,
                          input_file_name,
                          outfile,
                          MAX_SEGMENT_LENGTH);
            }
        }
        
        
        
        else if (strcmp(argv[1], "c_curve") == 0) {
            
            /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
            OptionParser opt_parse(strip_path(argv[1]),
                                   "", "<sorted-bed-file>");
            opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                              false , outfile);
            opt_parse.add_opt("step",'s',"step size in extrapolations "
                              "(default: " + toa(step_size) + ")",
                              false, step_size);
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
            opt_parse.add_opt("hist", 'H',
                              "input is a text file containing the observed histogram",
                              false, HIST_INPUT);
            opt_parse.add_opt("vals", 'V',
                              "input is a text file containing only the observed counts",
                              false, VALS_INPUT);
            
            
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
            /******************************************************************/
            
            if (argc > 2){
                // cerr << leftover_args.front();
                c_curve (VERBOSE,
                         VALS_INPUT,
                         PAIRED_END,
                         HIST_INPUT,
#ifdef HAVE_SAMTOOLS
                         BAM_FORMAT_INPUT,
#endif
                         step_size,
                         upper_limit,
                         input_file_name,
                         outfile,
                         MAX_SEGMENT_LENGTH);
            }
        }
        
        
        
        
        /*
         else if (strcmp(argv[1], "gc_extrap") == 0) {
         
         // ********* GET COMMAND LINE ARGUMENTS  FOR GC EXTRAP **********
         OptionParser opt_parse(strip_path(argv[1]),
         "", "<sorted-mapped-read-file>");
         opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
         false , outfile);
         opt_parse.add_opt("max_width", 'w', "max fragment length, "
         "set equal to read length for single end reads",
         false, max_width);
         opt_parse.add_opt("bin_size", 'n', "bin size "
         "(default: " + toa(bin_size) + ")",
         false, bin_size);
         opt_parse.add_opt("extrap",'e',"maximum extrapolation in base pairs"
         "(default: " + toa(max_extrapolation) + ")",
         false, max_extrapolation);
         //changed this to step size, originally was read step size. 
         opt_parse.add_opt("step",'s',"step size in bases between extrapolations "
         "(default: " + toa(step_size) + ")",
         false, step_size);
         opt_parse.add_opt("bootstraps",'b',"number of bootstraps "
         "(default: " + toa(bootstraps) + "), ",
         false, bootstraps);
         opt_parse.add_opt("cval", 'c', "level for confidence intervals "
         "(default: " + toa(c_level) + ")", false, c_level);
         opt_parse.add_opt("reads_per_base", 'r', "average reads per base "
         "(including duplicates) to predict sequencing level",
         false, reads_per_base);
         opt_parse.add_opt("terms",'x',"maximum number of terms",
         false, orig_max_terms);
         opt_parse.add_opt("fixed_fold",'f',"fixed fold extrapolation to predict",
         false, fixed_fold);
         //    opt_parse.add_opt("tol",'t', "numerical tolerance", false, tolerance);
         //    opt_parse.add_opt("max_iter",'i', "maximum number of iteration",
         //		      false, max_iter);
         opt_parse.add_opt("verbose", 'v', "print more information", 
         false, VERBOSE);
         opt_parse.add_opt("bed", 'B', "input is in bed format without sequence information",
         false, NO_SEQUENCE);
         opt_parse.add_opt("quick",'Q', 
         "quick mode: run gc_extrap without bootstrapping for confidence intervals",
         false, SINGLE_ESTIMATE);
         
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
         
         if (argc > 2){
         gc_extrap(VERBOSE,
         NO_SEQUENCE,
         SINGLE_ESTIMATE,
         MIN_REQUIRED_COUNTS,
         orig_max_terms,
         max_extrapolation,
         step_size,
         bootstraps,
         diagonal,
         c_level,
         max_width,
         bin_size,
         max_iter,
         tolerance,
         reads_per_base,
         fixed_fold,
         input_file_name,
         outfile);
         }
         
         }
         */
        
        else {
            cerr << "unrecognized command " << argv[1] << endl;
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
