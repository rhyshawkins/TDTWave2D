#pragma once
#ifndef value_hpp
#define value_hpp

#include <mpi.h>

#include "global.hpp"

class Value {
public:

  Value(Global &global);
  ~Value();

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

  int choose_value_location_and_value(int &value_depth,
				      int &value_idx,
				      double &choose_prob,
				      double &value,
				      int &ii,
				      int &ij,
				      double &value_prior_ratio,
				      int &prior_errors,
				      int &valid_proposal);

  int communicate_value_location_and_value(int &valid_proposal,
					   int &value_idx,
					   int &value_depth,
					   double &value);

  int propose_value(int valid_proposal,
		    int value_idx,
		    int value_depth,
		    double value);

  int compute_likelihood(int valid_idx, double &proposed_likelihood, double &proposed_log_normalization);

  int compute_acceptance(int value_idx,
			 double value_prior_ratio,
			 double proposed_likelihood,
			 double proposed_log_normalization,
			 bool &accept_proposal);
  
  int communicate_acceptance(bool &accept_proposal);
};

#endif // value_hpp
