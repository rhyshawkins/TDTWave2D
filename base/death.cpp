
#include <math.h>

extern "C" {
#include "slog.h"
};

#include "death.hpp"

#include "tdtwave2dexception.hpp"
#include "tdtwave2dutil.hpp"

Death::Death(Global &_global) :
  global(_global),
  propose(0),
  accept(0),
  propose_depth(new int[global.treemaxdepth + 1]),
  accept_depth(new int[global.treemaxdepth + 1]),
  communicator(MPI_COMM_NULL),
  mpi_size(-1),
  mpi_rank(-1)
{
  for (int i = 0; i <= global.treemaxdepth; i ++) {
    propose_depth[i] = 0;
    accept_depth[i] = 0;
  }
}

Death::~Death()
{
  delete [] propose_depth;
  delete [] accept_depth;
}

int
Death::step()
{
  propose ++;

  int k = wavetree2d_sub_coeff_count(global.wt);

  if (wavetree2d_sub_set_invalid_perturbation(global.wt, WT_PERTURB_DEATH) < 0) {
    return -1;
  }

  if (k > 1) {

    double ratio;
    int death_valid = 1;
    int death_depth;
    int death_idx;
    double choose_prob;
    double death_value;
    int ii, ij;
    double death_parent_coeff;
    double death_prob;
    double reverse_prob;
    double prior_prob;
    double proposed_likelihood;
    double proposed_log_normalization;
    
    //
    // First determine coefficient to death
    //
    if (choose_death_location(k,
			      ratio,
			      death_depth,
			      death_idx,
			      choose_prob,
			      death_valid) < 0) {
      return -1;
    }

    if (communicate_death_location(death_valid,
				   death_idx,
				   death_depth) < 0) {
      return -1;
    }

    
    if (death_valid) {

      propose_depth[death_depth] ++;
      
      if (propose_death(death_valid,
			death_idx,
			death_depth,
			death_value) < 0) {
	return -1;
      }
      
      if (compute_reverse_death_probability(death_idx,
					    death_depth,
					    death_value,
					    ii,
					    ij,
					    death_parent_coeff,
					    death_prob,
					    reverse_prob,
					    prior_prob) < 0) {
	return -1;
      }

      if (compute_likelihood(death_idx,
			     proposed_likelihood,
			     proposed_log_normalization) < 0) {
	return -1;
      }

      bool accept_proposal = false;

      if (compute_acceptance(proposed_likelihood,
			     proposed_log_normalization,
			     reverse_prob,
			     choose_prob,
			     death_prob,
			     ratio,
			     prior_prob,
			     accept_proposal) < 0) {
	return -1;
      }

      if (communicate_acceptance(accept_proposal) < 0) {
	return -1;
      }
      
      /*
       * Accept
       */
      if (accept_proposal) {
	accept ++;
	accept_depth[death_depth] ++;
	
        if (coefficient_histogram_accept_death(global.coeff_hist, death_idx) < 0) {
          ERROR("failed to update histogram for death acceptance\n");
          return -1;
        }
	
        if (wavetree2d_sub_commit(global.wt) < 0) {
          ERROR("failed to commit death\n");
          return -1;
        }
        
        global.current_likelihood = proposed_likelihood;
	global.current_log_normalization = proposed_log_normalization;
	global.accept();

        return 1;

      } else {
        
        /*
         * Reject
         */
        if (wavetree2d_sub_undo(global.wt) < 0) {
          ERROR("failed to undo death\n");
          return -1;
        }

	global.reject();
	
        return 0;
      }
    }
  }
	
  return 0;
}

std::string
Death::write_short_stats()
{
  return mkformatstring("Death %6d/%6d %7.3f",
			accept,
			propose,
			propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
}

std::string
Death::write_long_stats()
{
  std::string s = mkformatstring("Death: %6d %7.3f:",
				 propose,
				 propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
  
  for (int i = 0; i <= global.treemaxdepth; i ++) {
    s = s + mkformatstring("%7.3f ",
			   propose_depth[i] == 0 ? 0.0 : 100.0*(double)accept_depth[i]/(double)propose_depth[i]);
  }

  return s;
}

void
Death::initialize_mpi(MPI_Comm _communicator)
{
  MPI_Comm_dup(_communicator, &communicator);

  if (MPI_Comm_size(communicator, &mpi_size) != MPI_SUCCESS) {
    throw TDTWAVE2DEXCEPTION("MPI Failure\n");
  }
  if (MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS) {
    throw TDTWAVE2DEXCEPTION("MPI Failure\n");
  }
}

bool
Death::primary() const
{
  return (communicator == MPI_COMM_NULL || mpi_rank == 0);
}

int
Death::choose_death_location(int k,
			     double &ratio,
			     int &death_depth,
			     int &death_idx,
			     double &choose_prob,
			     int &death_valid)
{
  if (primary()) {
    //
    // Get tree structure ratio
    //
    
    if (hnk_get_kplus1_ratio(global.hnk, 
			     global.treemaxdepth, 
			     k - 1, 
			     &ratio) < 0) {
      ERROR("failed to get ratio for birth\n");
      return -1;
    }
    
    if (wavetree2d_sub_choose_death_global(global.wt, 
					   global.random.uniform(), 
					   global.treemaxdepth, 
					   &death_depth,
					   &death_idx,
					   &choose_prob) < 0) {
      death_valid = 0;
    }
  }
  
  return 0;
}

int
Death::communicate_death_location(int &death_valid,
				  int &death_idx,
				  int &death_depth)
{
  if (communicator != MPI_COMM_NULL) {

    if (MPI_Bcast(&death_valid, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw TDTWAVE2DEXCEPTION("Failed to broadcast death valid\n");
    }
    
    if (death_valid) {
      if (MPI_Bcast(&death_idx, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
	throw TDTWAVE2DEXCEPTION("Failed to broadcast death index\n");
      }
      if (MPI_Bcast(&death_depth, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
	throw TDTWAVE2DEXCEPTION("Failed to broadcast death depth\n");
      }
    }
  }

  return 0;
}

int
Death::propose_death(int death_valid,
		     int death_idx,
		     int death_depth,
		     double &death_value)
{
  if (death_valid) {
    if (coefficient_histogram_propose_death(global.coeff_hist, death_idx) < 0) {
      ERROR("failed to update histogram for death proposal\n");
      return -1;
    }
    
    if (wavetree2d_sub_propose_death(global.wt,
				     death_idx, 
				     death_depth,
				     &death_value) < 0) {
      ERROR("failed to propose death\n");
      return -1;
    }
  }

  return 0;
}

int
Death::compute_reverse_death_probability(int death_idx,
					 int death_depth,
					 double death_value,
					 int &ii,
					 int &ij,
					 double &death_parent_coeff,
					 double &death_prob,
					 double &reverse_prob,
					 double &prior_prob)
{
  if (primary()) {
    
    if (wavetree2d_sub_2dindices(global.wt, death_idx, &ii, &ij) < 0) {
      ERROR("failed to compute 2d indices for death\n");
      return -1;
    }
    
    if (wavetree2d_sub_get_coeff(global.wt,
				 wavetree2d_sub_parent_index(global.wt, death_idx),
				 &death_parent_coeff) < 0) {
      ERROR("failed to get parent coefficient for death\n");
      return -1;
    }
    
    if (wavetree_pp_death2d(global.proposal,
			    ii,
			    ij,
			    death_depth,
			    global.treemaxdepth,
			    death_parent_coeff,
			    death_value,
			    &death_prob) < 0) {
      ERROR("failed to get death probability\n");
      return -1;
    }
    
    if (wavetree2d_sub_reverse_death_global(global.wt, 
					    global.treemaxdepth,
					    death_depth,
					    death_idx, 
					    &reverse_prob) < 0) {
      ERROR("failed to reverse death (global)\n");
      return -1;
    }
    
    //
    // Compute the prior ratio
    //
    prior_prob = wavetree_pp_prior_probability2d(global.proposal,
						 ii, 
						 ij,
						 death_depth,
						 global.treemaxdepth,
						 death_parent_coeff,
						 death_value);
  }

  return 0;
}

int
Death::compute_likelihood(int death_idx,
			  double &proposed_likelihood,
			  double &proposed_log_normalization)
{
  if (communicator == MPI_COMM_NULL) {

    proposed_likelihood = global.likelihood(proposed_log_normalization);

  } else {

    proposed_likelihood = global.likelihood_mpi(proposed_log_normalization);

  }

  return 0;
}

int
Death::compute_acceptance(double proposed_likelihood,
			  double proposed_log_normalization,
			  double reverse_prob,
			  double choose_prob,
			  double death_prob,
			  double ratio,
			  double prior_prob,
			  bool &accept_proposal)
{
  if (primary()) {

    double u = log(global.random.uniform());
    
    accept_proposal = u < ((global.current_likelihood - proposed_likelihood)/global.temperature // Likelihood ratio
			   + log(reverse_prob) - 
			   log(choose_prob)                           // Node proposal ratio 
			   + log(death_prob)                          // Coefficient proposal ratio 
			   + log(ratio)                               // Tree prior ratio 
			   - log(prior_prob)                          // Coefficient prior ratio 
			   );

  }

  return 0;
}

int
Death::communicate_acceptance(bool &accept_proposal)
{
  if (communicator != MPI_COMM_NULL) {

    int ta;

    if (mpi_rank == 0) {
      ta = (int)accept_proposal;
    }

    if (MPI_Bcast(&ta, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw TDTWAVE2DEXCEPTION("Failed to broadcast acceptance\n");
    }

    if (mpi_rank != 0) {
      accept_proposal = (bool)ta;
    }
  }

  return 0;
}
