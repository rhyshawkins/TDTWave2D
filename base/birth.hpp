#pragma once
#ifndef birth_hpp
#define birth_hpp

#include "global.hpp"

#include <mpi.h>

class Birth {
public:

  Birth(Global &global);
  ~Birth();

  int step();

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm communicator);

  Global &global;
  int propose;
  int accept;

  int *propose_depth;
  int *accept_depth;

  MPI_Comm communicator;
  int mpi_size;
  int mpi_rank;

private:

  bool primary() const;

  int choose_birth_location_and_value(int k,
				      double &ratio,
				      int &birth_depth,
				      int &birth_idx,
				      double &choose_prob,
				      double &birth_value,
				      double &birth_prob,
				      int &birth_valid,
				      double &birth_parent_coeff,
				      int &ii,
				      int &ij);

  int communicate_birth_location_and_value(int &birth_valid,
					   int &birth_idx,
					   int &birth_depth,
					   double &birth_value);

  int propose_birth(int &birth_valid,
		    int &birth_idx,
		    int &birth_depth,
		    double &birth_value);

  int compute_reverse_birth_probability(int birth_depth,
					int birth_idx,
					int ii,
					int ij,
					double birth_parent_coeff,
					double birth_value,
					double &reverse_prob,
					double &prior_prob);

  int compute_likelihood(int birth_idx, double &proposed_likelihood, double &proposed_log_normalization);

  int compute_acceptance(double proposed_likelihood,
			 double proposed_log_normalization,
			 double reverse_prob,
			 double choose_prob,
			 double birth_prob,
			 double ratio,
			 double prior_prob,
			 bool &accept_proposal);

  int communicate_acceptance(bool &accept_proposal);

};

#endif // birth_hpp
