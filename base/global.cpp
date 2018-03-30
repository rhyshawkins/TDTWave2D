
#include "tdtwave2dexception.hpp"

#include "global.hpp"

extern "C" {
  #include "hnk_cartesian_nonsquare.h"

  #include "slog.h"
};

#include "constants.hpp"

const int CHAIN_STEPS = 1000000;

int global_coordtoindex(void *user, int i, int j, int k, int depth)
{
  wavetree2d_sub_t *wt = (wavetree2d_sub_t*)user;

  return wavetree2d_sub_from_2dindices(wt, i, j);
}

int global_indextocoord(void *user, int index, int *i, int *j, int *k, int *depth)
{
  wavetree2d_sub_t *wt = (wavetree2d_sub_t*)user;

  if (wavetree2d_sub_2dindices(wt, index, i, j) < 0) {
    return -1;
  }

  *k = 0;
  *depth = wavetree2d_sub_depthofindex(wt, index);
  return 0;
}

Global::Global(const char *filename,
	       const char *initial_model,
	       const char *prior_file,
	       const char *hierarchical_prior_file,
	       int _degreex,
	       int _degreey,
	       int seed,
	       int _kmax,
	       bool _posteriork,
	       int hwavelet,
	       int vwavelet) :
  kmax(_kmax),
  treemaxdepth(-1),
  wt(nullptr),
  ch(nullptr),
  hnk(nullptr),
  proposal(nullptr),
  degreex(_degreex),
  degreey(_degreey),
  image(nullptr),
  model(nullptr),
  workspace(nullptr),
  mean_residual_n(0),
  residual(nullptr),
  mean_residual(nullptr),
  last_valid_residual(nullptr),
  residual_normed(nullptr),
  mean_residual_normed(nullptr),
  last_valid_residual_normed(nullptr),
  residuals_valid(false),
  residual_hist_bins(100),
  residual_hist_min(-5.0),
  residual_hist_max(5.0),
  residual_hist(nullptr),
  width(-1),
  height(-1),
  size(-1),
  ncoeff(-1),
  lambda_scale(1.0),
  current_likelihood(-1.0),
  coeff_hist(nullptr),
  random(seed),
  posteriork(_posteriork),
  communicator(MPI_COMM_NULL),
  mpi_size(-1),
  mpi_rank(-1),
  temperature(1.0),
  column_offsets(nullptr),
  column_sizes(nullptr),
  residual_offsets(nullptr),
  residual_sizes(nullptr)
{
  if (degreex < 0 || degreex >= 16 ||
      degreey < 0 || degreey >= 16) {
    throw TDTWAVE2DEXCEPTION("Degree(s) out of range: %d x %d\n", degreex, degreey);
  }
  
  if (depth <= 0.0) {
    throw TDTWAVE2DEXCEPTION("Depth out of range\n");
  }

  if (!posteriork) {
    //
    // Load observations
    //

  }

  wt = wavetree2d_sub_create(degreex, degreey, 0.0);
  if (wt == NULL) {
    throw TDTWAVE2DEXCEPTION("Failed to create wavetree\n");
  }

  width = wavetree2d_sub_get_width(wt);
  height = wavetree2d_sub_get_height(wt);
  size = wavetree2d_sub_get_size(wt);
  ncoeff = wavetree2d_sub_get_ncoeff(wt);
  treemaxdepth = wavetree2d_sub_maxdepth(wt);

  INFO("Image: %d x %d\n", width, height);

  if (!posteriork) {
    
    image = new tdtwave2dimage(height, width, 0.0);

    model = new double[size];
    int workspacesize = width;
    if (height > width) {
      workspacesize = height;
    }
    workspace = new double[workspacesize];

    int ntotal = 0; // TODO
    INFO("Data: %d total points\n", ntotal);

    residual_size = ntotal;
    residual = new double[residual_size];
    mean_residual = new double[residual_size];
    last_valid_residual = new double[residual_size];
    residual_normed = new double[residual_size];
    mean_residual_normed = new double[residual_size];
    last_valid_residual_normed = new double [residual_size];

    residuals_per_column = residual_size/image->columns;

    residual_hist = new int[residual_size * residual_hist_bins];

    reset_residuals();
  }
  
  if (initial_model == NULL) {
    if (wavetree2d_sub_initialize(wt, 0.0) < 0) {
      throw TDTWAVE2DEXCEPTION("Failed to initialize wavetree\n");
    }
  } else {

    if (wavetree2d_sub_load_promote(wt, initial_model) < 0) {
      throw TDTWAVE2DEXCEPTION("Failed to load initial model: %s\n", initial_model);
    }

    INFO("Loaded model with %d coefficients\n", wavetree2d_sub_coeff_count(wt));
  }

  //
  // Hnk Ratio
  //
  if (kmax > ncoeff) {
    INFO("Warning: kmax truncated to %d\n", ncoeff);
    kmax = ncoeff;
  }
  
  hnk = hnk_cartesian_nonsquare_2D_create_sub(degreex,
					      degreey,
					      kmax);
  if (hnk == NULL) {
    throw TDTWAVE2DEXCEPTION("Failed to create hnk table\n");
  }

  //
  // Chain History
  //
  ch = chain_history_create(CHAIN_STEPS);
  if (ch == nullptr) {
    throw TDTWAVE2DEXCEPTION("Failed to create chain history\n");
  }

  
  //
  // Initialse coeff histogram
  //
  coeff_hist = coefficient_histogram_create(ncoeff, 100, -1.0, 1.0,
					    global_coordtoindex,
					    global_indextocoord,
					    wt);
  if (coeff_hist == NULL) {
    throw TDTWAVE2DEXCEPTION("Failed to create coefficient histogram\n");
  }

  //
  // Create proposal structure.
  //
  if (prior_file != nullptr) {
    proposal = wavetree_pp_load(prior_file, seed, coeff_hist);
    if (proposal == NULL) {
      throw TDTWAVE2DEXCEPTION("Failed to load proposal file\n");
    }
    
    for (int i = 0; i < ncoeff; i ++) {
      int depth = wavetree2d_sub_depthofindex(wt, i);
      
      int ii, ij;
      if (wavetree2d_sub_2dindices(wt, i, &ii, &ij) < 0) {
	throw TDTWAVE2DEXCEPTION("Failed to get 2d indices\n");
      }
      
      double vmin, vmax;
      if (wavetree_pp_prior_range2d(proposal,
				    ii,
				    ij,
				    depth,
				    treemaxdepth,
				    0.0,
				    &vmin,
				    &vmax) < 0) {
	throw TDTWAVE2DEXCEPTION("Failed to get coefficient range\n");
      }
      
      if (coefficient_histogram_set_range(coeff_hist, 
					  i,
					  vmin,
					  vmax) < 0) {
	throw TDTWAVE2DEXCEPTION("Failed to set coefficient histogram range\n");
      }
    }
  }
    
  hwaveletf = wavelet_inverse_function_from_id(hwavelet);
  if (hwaveletf == nullptr) {
    throw TDTWAVE2DEXCEPTION("Invalid horizontal wavelet %d\n", hwavelet);
  }

  vwaveletf = wavelet_inverse_function_from_id(vwavelet);
  if (vwaveletf == nullptr) {
    throw TDTWAVE2DEXCEPTION("Invalid vertical wavelet %d\n", vwavelet);
  }

}

Global::~Global()
{
}

double
Global::likelihood(double &log_normalization)
{
  if (!posteriork) {
    
    //
    // Get tree model wavelet coefficients
    //
    memset(image->image, 0, sizeof(double) * size);
    if (wavetree2d_sub_map_to_array(wt, image->image, size) < 0) {
      throw TDTWAVE2DEXCEPTION("Failed to map model to array\n");
    }

    //
    // Inverse wavelet transform
    //
    if (generic_lift_inverse2d(image->image,
			       width,
			       height,
			       width,
			       workspace,
			       hwaveletf,
			       vwaveletf,
			       1) < 0) {
      throw TDTWAVE2DEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    double sum = 0.0;
    int residual_offset;

    log_normalization = 0.0;
    
    for (int i = 0; i < image->columns; i ++) {

      residual_offset = i * residuals_per_column;
      

      
    }
    
    return sum;
  } else {
    return 1.0;
  }
}

double
Global::hierarchical_likelihood(double proposed_lambda_scale,
				double &log_normalization)
{
  log_normalization = 0.0;

  if (!posteriork) {

    if (!residuals_valid) {
      double x;
      (void)likelihood(x);
      accept();
    }

    double sum = 0.0;
    int residual_offset;
  
    for (int i = 0; i < image->columns; i ++) {
      
      residual_offset = i * residuals_per_column;
    }

    return sum;

  } else {
    return 1.0;
  }
}

void
Global::initialize_mpi(MPI_Comm _communicator, double _temperature)
{
  communicator = _communicator;

  if (MPI_Comm_size(communicator, &mpi_size) != MPI_SUCCESS) {
    throw TDTWAVE2DEXCEPTION("MPI Failure\n");
  }
  
  if (MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS) {
    throw TDTWAVE2DEXCEPTION("MPI Failure\n");
  }

  column_offsets = new int[mpi_size];
  column_sizes = new int[mpi_size];
  residual_offsets = new int[mpi_size];
  residual_sizes = new int[mpi_size];

  int columns = image->columns;
  int processes = mpi_size;

  //
  // Evenly distribute columns
  //
  for (int i = 0; i < mpi_size; i ++) {
    column_sizes[i] = columns/processes;
    residual_sizes[i] = column_sizes[i] * residuals_per_column;

    columns -= column_sizes[i];
    processes --;
  }

  column_offsets[0] = 0;
  residual_offsets[0] = 0;
  for (int i = 1; i < mpi_size; i ++) {
    column_offsets[i] = column_offsets[i - 1] + column_sizes[i - 1];
    residual_offsets[i] = column_offsets[i] * residuals_per_column;
    INFO("Split: %4d %4d", column_offsets[i], column_sizes[i]);
  }

  if (column_offsets[mpi_size - 1] + column_sizes[mpi_size - 1] != image->columns) {
    throw TDTWAVE2DEXCEPTION("Column sharing intialization failure\n");
  }

  temperature = _temperature;
}

double
Global::likelihood_mpi(double &log_normalization)
{
  if (communicator == MPI_COMM_NULL || mpi_rank < 0 || mpi_size < 0) {
    throw TDTWAVE2DEXCEPTION("MPI Parameters unset\n");
  }
  
  if (!posteriork) {
    
    //
    // Get tree model wavelet coefficients
    //
    memset(image->image, 0, sizeof(double) * size);
    if (wavetree2d_sub_map_to_array(wt, image->image, size) < 0) {
      throw TDTWAVE2DEXCEPTION("Failed to map model to array\n");
    }

    //
    // Inverse wavelet transform
    //
    if (generic_lift_inverse2d(image->image,
			       width,
			       height,
			       width,
			       workspace,
			       hwaveletf,
			       vwaveletf,
			       1) < 0) {
      throw TDTWAVE2DEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    double sum = 0.0;
    double local_log_normalization = 0.0;
    int residual_offset;
    
    for (int mi = 0, i = column_offsets[mpi_rank]; mi < column_sizes[mpi_rank]; mi ++, i ++) {

      residual_offset = i * residuals_per_column;
      
    }

    double total;
    
    if (MPI_Reduce(&local_log_normalization, &total, 1, MPI_DOUBLE, MPI_SUM, 0, communicator) != MPI_SUCCESS) {
      throw TDTWAVE2DEXCEPTION("Likelihood failed in reducing\n");
    }
    if (MPI_Bcast(&total, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
      throw TDTWAVE2DEXCEPTION("Likelihood failed in broadcast\n");
    }

    log_normalization = total;
    
    if (MPI_Reduce(&sum, &total, 1, MPI_DOUBLE, MPI_SUM, 0, communicator) != MPI_SUCCESS) {
      throw TDTWAVE2DEXCEPTION("Likelihood failed in reducing\n");
    }
    if (MPI_Bcast(&total, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
      throw TDTWAVE2DEXCEPTION("Likelihood failed in broadcast\n");
    }


    MPI_Allgatherv(residual + residual_offsets[mpi_rank],
		   residual_sizes[mpi_rank],
		   MPI_DOUBLE,
		   residual,
		   residual_sizes,
		   residual_offsets,
		   MPI_DOUBLE,
		   communicator);

    MPI_Allgatherv(residual_normed + residual_offsets[mpi_rank],
		   residual_sizes[mpi_rank],
		   MPI_DOUBLE,
		   residual_normed,
		   residual_sizes,
		   residual_offsets,
		   MPI_DOUBLE,
		   communicator);

    return total;
    
  } else {
    log_normalization = 0.0;
    
    return 1.0;
  }
}

double
Global::hierarchical_likelihood_mpi(double proposed_lambda_scale,
				    double &log_normalization)
{
  double like;

  log_normalization = 0.0;

  if (!residuals_valid) {
    double x;
    like = likelihood_mpi(x);
    accept();
  }
  
  //
  // 
  //
  // Since this is just a sum over residuals, just compute on root node and broadcast.
  //
  if (mpi_rank == 0) {
    like = hierarchical_likelihood(proposed_lambda_scale, log_normalization);
  }
  
  MPI_Bcast(&like, 1, MPI_DOUBLE, 0, communicator);
  MPI_Bcast(&log_normalization, 1, MPI_DOUBLE, 0, communicator);

  return like;
}

void
Global::reset_residuals()
{
  mean_residual_n = 0;
  for (int i = 0; i < residual_size; i ++) {
    residual[i] = 0.0;
    mean_residual[i] = 0.0;
    last_valid_residual[i] = 0.0;
      
    residual_normed[i] = 0.0;
    mean_residual_normed[i] = 0.0;
    last_valid_residual_normed[i] = 0.0;

    for (int j = 0; j < residual_hist_bins; j ++) {
      residual_hist[i * residual_hist_bins + j] = 0;
    }
  }

  cov_n = 0;

  for (int i = 0; i < (int)cov_count.size(); i ++) {
    int N = cov_count[i];

    for (int j = 0; j < N; j ++) {
      cov_delta[i][j] = 0.0;
      cov_mu[i][j] = 0.0;
    }

    N = N*N;
    for (int j = 0; j < N; j ++) {
      cov_sigma[i][j] = 0.0;
    }
  }
}

void
Global::invalidate_residuals()
{
  residuals_valid = false;
}

void
Global::accept()
{
  residuals_valid = true;
  if (!posteriork) {
    for (int i = 0; i < residual_size; i ++) {
      last_valid_residual[i] = residual[i];
      last_valid_residual_normed[i] = residual_normed[i];
    }

    update_residual_mean();
    update_residual_covariance();
  }
}

void
Global::accept_hierarchical()
{
}

void
Global::reject()
{
  update_residual_mean();
}

void
Global::reject_hierarchical()
{
}

void
Global::update_residual_mean()
{
  mean_residual_n ++;

  for (int i = 0; i < residual_size; i ++) {

    double delta = last_valid_residual[i] - mean_residual[i];
    mean_residual[i] += delta/(double)(mean_residual_n);

    delta = last_valid_residual_normed[i] - mean_residual_normed[i];
    mean_residual_normed[i] += delta/(double)(mean_residual_n);

    int hi = (int)((last_valid_residual_normed[i] - residual_hist_min)/(residual_hist_max - residual_hist_min) * (double)residual_hist_bins);
    if (hi >= 0 && hi < residual_hist_bins) {
      residual_hist[i * residual_hist_bins + hi] ++;
    }
  }
}

void
Global::update_residual_covariance()
{
  double *p = last_valid_residual;
  int nobservations = 0; //TODO
  for (int k = 0; k < nobservations; k ++) {

    cov_n ++;
    
    for (int i = 0; i < (int)cov_count.size(); i ++) {

      int N = cov_count[i];
      
      for (int j = 0; j < N; j ++) {
	cov_delta[i][j] = (p[j] - cov_mu[i][j])/(double)(cov_n);
	cov_mu[i][j] += cov_delta[i][j];
      }

      for (int j = 0; j < N; j ++) {
	for (int l = j; l < N; l ++) {

	  cov_sigma[i][j * N + l] +=
	    (double)(cov_n - 1)*cov_delta[i][j]*cov_delta[i][l] -
	    cov_sigma[i][j * N + l]/(double)(cov_n);
	}
      }

      p += N;
    }
  }
}


int
Global::get_residual_size() const
{
  return residual_size;
}

const double *
Global::get_mean_residuals() const
{
  return mean_residual;
}
  
const double *
Global::get_mean_normed_residuals() const
{
  return mean_residual_normed;
}

bool
Global::save_residual_histogram(const char *filename) const
{
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("Failed to create file");
    return false;
  }

  fprintf(fp, "%d %d %f %f\n", residual_size, residual_hist_bins, residual_hist_min, residual_hist_max);
  for (int i = 0; i < residual_size; i ++) {
    for (int j = 0; j < residual_hist_bins; j ++) {
      fprintf(fp, "%d ", residual_hist[i * residual_hist_bins + j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  return true;
}

bool
Global::save_residual_covariance(const char *filename) const
{
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("Failed to create file\n");
    return false;
  }

  fprintf(fp, "%d\n", (int)cov_count.size());
  
  for (int i = 0; i < (int)cov_count.size(); i ++) {

    int N = cov_count[i];

    fprintf(fp, "%d\n", N);
    for (int j = 0; j < N; j ++) {
      fprintf(fp, "%.9g ", cov_mu[i][j]);
    }
    fprintf(fp, "\n");

    for (int j = 0; j < N; j ++) {
      for (int k = 0; k < N; k ++) {
	fprintf(fp, "%.9g ", cov_sigma[i][j * N + k]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);

  return true;
}

generic_lift_inverse1d_step_t
Global::wavelet_inverse_function_from_id(int id)
{
  switch (id) {
  case WAVELET_HAAR:
    return haar_lift_inverse1d_haar_step;

  case WAVELET_DAUB4:
    return daub4_dwt_inverse1d_daub4_step;

  case WAVELET_DAUB6:
    return daub6_dwt_inverse1d_daub6_step;

  case WAVELET_DAUB8:
    return daub8_dwt_inverse1d_daub8_step;

  case WAVELET_CDF97:
    return cdf97_lift_inverse1d_cdf97_step;

  case WAVELET_CDF97_PERIODIC:
    return cdf97_lift_periodic_inverse1d_cdf97_step;

  default:
    return nullptr;
  }
}

generic_lift_forward1d_step_t
Global::wavelet_forward_function_from_id(int id)
{
  switch (id) {
  case WAVELET_HAAR:
    return haar_lift_forward1d_haar_step;

  case WAVELET_DAUB4:
    return daub4_dwt_forward1d_daub4_step;

  case WAVELET_DAUB6:
    return daub6_dwt_forward1d_daub6_step;

  case WAVELET_DAUB8:
    return daub8_dwt_forward1d_daub8_step;

  case WAVELET_CDF97:
    return cdf97_lift_forward1d_cdf97_step;

  case WAVELET_CDF97_PERIODIC:
    return cdf97_lift_periodic_forward1d_cdf97_step;

  default:
    return nullptr;
  }
}
