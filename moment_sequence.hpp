/*    Copyright (C) 2013 University of Southern California and
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

#ifndef MOMENT_SEQUENCE_HPP
#define MOMENT_SEQUENCE_HPP

#include <vector>
#include <numeric>

// test Hankel moment matrix to ensure the moment sequence
// is positive definite
size_t ensure_pos_def_mom_seq(std::vector<double> &moments,
			      const double tolerance,
			      const bool VERBOSE);

struct MomentSequence {

  // Constructors
  MomentSequence() {}
  MomentSequence(const std::vector<double> &obs_moms);

  MomentSequence(const std::vector<double> &a,
		 const std::vector<double> &b):
    alpha(a), beta(b) {};



  // Estimate 3-term recurrence
  // these will be removed from the header when they are tested
  void gw_three_term_calc(const bool VERBOSE, const size_t n_points);
  void unmodified_Chebyshev(const bool VERBOSE);
  void modified_Chebyshev(const bool VERBOSE,
			  const size_t n_points,
			  const std::vector<double> &mod_alpha,
			  const std::vector<double> &mod_beta,
			  const std::vector<double> &modified_moments);

  void full_3term_recurrence(const bool VERBOSE,
			     std::vector<double> &full_alpha,
			     std::vector<double> &full_beta);

  // quadrature rules using polynomial solver
  void poly_solve_gauss_quad(const size_t n_points,
			     std::vector<double> &weights,
			     std::vector<double> &points);

  // quadrature rules using QR on Jacobi matrix 
  bool Lower_quadrature_rules(const bool VERBOSE,
				 const size_t n_points,
				 const double tolerance, 
				 const size_t max_iter,
				 std::vector<double> &points,
				 std::vector<double> &weights);

  bool GaussRadau_quadrature_rules(const bool VERBOSE,
				   const size_t n_points,
				   const double tolerance,
				   const size_t max_iter,
				   const double fixed_left_end_point,
				   std::vector<double> &points,
				   std::vector<double> &weights);

  // points are determined assuming data is NegBin distrbuted
  // 3-term recurrence is therefore known
  // weights are determined by satisfying observed moment conditions
  bool NegBin_quadrature_rules(const bool VERBOSE,
			       const size_t n_points,
			       const double tolerance,
			       const size_t max_iter,
			       const double estimated_mu,
			       const double estimated_alpha,
			       std::vector<double> &points,
			       std::vector<double> &weights);


  std::vector<double> moments;
  // 3-term recurrence
  std::vector<double> alpha;
  std::vector<double> beta;
};


#endif
