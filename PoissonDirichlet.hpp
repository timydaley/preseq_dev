/*    PoissonDirichlet:
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *                       Timothy Daley
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


#ifndef POISSONDIRICHLET_HPP
#define POISSONDIRICHLET_HPP

#include <fstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numeric>

#include <gsl/gsl_rng.h>


class PDD{
public:
  PDD(const double t_, const double s_):
    theta(t_), sigma(s_) {};
  double get_theta() const {return theta;}
  double get_sigma() const {return sigma;}
  
  void set_theta(const double t_) {theta = t_;}
  void set_sigma(const double s_) {sigma = s_;}

  void newton_raphson_estim_params(const bool VERBOSE,
				   const std::vector<double> &counts);

  double expected_additional_distinct(const size_t current_total_count,
				      const size_t current_distinct,
				      const size_t future_total_count);

  double expected_discovery_prob(const size_t current_total_count,
				 const size_t current_distinct,
				 const size_t future_total_count);

  void sample_PoissDir_counts(const gsl_rng *rng, const size_t total_count,
			      std::vector<double> &sampled_counts);

private:
  static const double min_allowed_sigma = 0.0;
  static const double max_allowed_sigma = 1.0;
  static const double min_allowed_theta = 0.0;
  static const double tolerance;
  static const size_t max_iter;

  double theta;
  double sigma;
};


#endif
