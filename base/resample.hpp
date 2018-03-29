#pragma once
#ifndef resample_hpp

#include <string>

#include "global.hpp"
#include <mpi.h>

class Resample {
public:

  Resample(Global &global);
  ~Resample();

  int step(double resample_temperature = 1.0);

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm global_communicator,
		      MPI_Comm temperature_communicator,
		      MPI_Comm chain_communicator);

  Global &global;

  MPI_Comm global_communicator;
  MPI_Comm temperature_communicator;
  MPI_Comm chain_communicator;

  int resamplings;
  int reselected;
  int propagated;
  int replaced;

  int global_size;
  int global_rank;

  int temperature_size;
  int temperature_rank;

  int chain_size;
  int chain_rank;

  int processesperchain;

  double *likelihood_array;
  double *lambda_array;
  double *weights_array;
  int *sources_array;
  MPI_Request *requests_array;

  int send_buffer_size;
  char *send_buffer;

  int recv_buffer_size;
  char *recv_buffer;

};

#endif // resample_hpp
