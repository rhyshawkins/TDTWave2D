#pragma once
#ifndef ptexchange_hpp
#define ptexchange_hpp

#include "global.hpp"
#include <mpi.h>

class PTExchange {
public:

  PTExchange(Global &global);
  ~PTExchange();

  int step();

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm global_communicator,
		      MPI_Comm temperature_communicator,
		      MPI_Comm chain_communicator,
		      int ntemperatures);

  Global &global;

  int propose;
  int accept;

  MPI_Comm global_communicator;
  int global_size;
  int global_rank;
  
  MPI_Comm temperature_communicator;
  int temperature_size;
  int temperature_rank;
  
  MPI_Comm chain_communicator;
  int chain_rank;

private:

  int ntotalchains;
  int processesperchain;
  int ntemperatures;
  int chainspertemperature;

  int *ptpairs;
  int *transposed_ptpairs;
  
  int partner;
  bool send;
  bool ptaccept;
  int send_length;
  int recv_length;
	
  double sendmsg[5];
  double recvmsg[5];
  
  double current_logprior;
  double partner_logprior;
  
  double partner_likelihood;
  double partner_log_normalization;
  double partner_temperature;
  double u;
  MPI_Status status;

  int exchanged;

  int send_buffer_size;
  char *send_buffer;

  int recv_buffer_size;
  char *recv_buffer;
  
  
};

#endif // ptexchange_hpp
