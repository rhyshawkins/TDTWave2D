#pragma once
#ifndef hierarchicalprior_hpp
#define hierarchicalprior_hpp

#include <mpi.h>

#include "global.hpp"

extern "C" {
  #include "chain_history.h"
};

class HierarchicalPrior {
public:

  HierarchicalPrior(Global &global, double sigma);
  ~HierarchicalPrior();

  int step();

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm communicator);

  void get_last_step(chain_history_change_t *last_step);

  Global &global;
  double sigma;

  int propose;
  int accept;

  MPI_Comm communicator;
  int mpi_size;
  int mpi_rank;

  chain_history_change_t last_step;

  double alpha;
  double beta;

private:

  bool primary() const;

  int choose_value(double &value,
		   double &value_prior_ratio,
		   int &valid_proposal);
  
  int communicate_value(int &valid_proposal,
			double &value);

  int compute_likelihood(double value,
			 double &proposed_likelihood,
			 double &proposed_log_normalization);

  int compute_acceptance(double log_value_prior_ratio,
			 double current_log_prior,
			 double proposed_log_prior,
			 bool &accept_proposal);
  
  int communicate_acceptance(bool &accept_proposal);

};

#endif // hierarchicalprior_hpp
