
#include "tdtwave2dexception.hpp"

#include "global.hpp"

extern "C" {
  #include "hnk_cartesian_nonsquare.h"

  #include "slog.h"
};

#include "constants.hpp"

#include "genericinterface.hpp"

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
	       int wavelet) :
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
  nobservations(0),
  predictions(nullptr),
  unused(nullptr),
  mean_residual_n(0),
  residual(nullptr),
  last_valid_residual(nullptr),
  mean_residual(nullptr),
  residual_normed(nullptr),
  last_valid_residual_normed(nullptr),
  mean_residual_normed(nullptr),
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
  residual_offsets(nullptr),
  residual_sizes(nullptr)
{
  if (degreex < 0 || degreex >= 16 ||
      degreey < 0 || degreey >= 16) {
    throw TDTWAVE2DEXCEPTION("Degree(s) out of range: %d x %d\n", degreex, degreey);
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

    int n = strlen(filename);
    nobservations = tdtwave2d_loaddata_(&n, filename, &width, &height);
    if (nobservations < 0) {
      throw TDTWAVE2DEXCEPTION("Failed to read input data file");
    }

    INFO("Data: %d total points\n", nobservations);

    residual_size = nobservations;
    predictions = new double[residual_size];
    unused = new double[residual_size];
    
    residual = new double[residual_size];
    last_valid_residual = new double[residual_size];
    mean_residual = new double[residual_size];
    residual_normed = new double[residual_size];
    last_valid_residual_normed = new double[residual_size];
    mean_residual_normed = new double[residual_size];

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
    
  hwaveletf = wavelet_inverse_function_from_id(wavelet);
  if (hwaveletf == nullptr) {
    throw TDTWAVE2DEXCEPTION("Invalid horizontal wavelet %d\n", wavelet);
  }

  vwaveletf = wavelet_inverse_function_from_id(wavelet);
  if (vwaveletf == nullptr) {
    throw TDTWAVE2DEXCEPTION("Invalid vertical wavelet %d\n", wavelet);
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

    for (int i = 0; i < nobservations; i ++) {
      if (tdtwave2d_compute_prediction_(&i,
					&width,
					&height,
					image->image,
					unused,
					&predictions[i]) < 0) {
	throw TDTWAVE2DEXCEPTION("Failed to compute prediction for observations %d\n", i);
      }
					
    }

    double likelihood = 0.0;
    double hvalue = lambda_scale;
    if (tdtwave2d_compute_likelihood_(&nobservations,
				      &hvalue,
				      predictions,
				      residual,
				      unused,
				      &likelihood,
				      &log_normalization) < 0) {
      throw TDTWAVE2DEXCEPTION("Failed to compute likelihood");
    }
    
    return likelihood;
    
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

    double likelihood = 0.0;
    if (tdtwave2d_compute_likelihood_(&nobservations,
				      &proposed_lambda_scale,
				      predictions,
				      residual,
				      unused,
				      &likelihood,
				      &log_normalization) < 0) {
      throw TDTWAVE2DEXCEPTION("Failed to compute likelihood");
    }
    
    return likelihood;
    
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

  residual_offsets = new int[mpi_size];
  residual_sizes = new int[mpi_size];

  int columns = nobservations;
  int processes = mpi_size;

  //
  // Evenly distribute columns
  //
  for (int i = 0; i < mpi_size; i ++) {
    residual_sizes[i] = columns/processes;

    columns -= residual_sizes[i];
    processes --;
  }

  residual_offsets[0] = 0;
  for (int i = 1; i < mpi_size; i ++) {
    residual_offsets[i] = residual_offsets[i - 1] + residual_sizes[i - 1];
    INFO("Split: %4d %4d", residual_offsets[i], residual_sizes[i]);
  }

  if (residual_offsets[mpi_size - 1] + residual_sizes[mpi_size - 1] != nobservations) {
    throw TDTWAVE2DEXCEPTION("Residual sharing intialization failure\n");
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

    //
    // In parallel, compute predictions
    //
    
    for (int mi = 0, i = residual_offsets[mpi_rank]; mi < residual_sizes[mpi_rank]; mi ++, i ++) {

      if (tdtwave2d_compute_prediction_(&i,
					&width,
					&height,
					image->image,
					unused,
					&predictions[i]) < 0) {
	throw TDTWAVE2DEXCEPTION("Failed to compute prediction for observations %d\n", i);
      }
					
    }

    //
    // Gather predictions
    //
    MPI_Allgatherv(predictions + residual_offsets[mpi_rank],
		   residual_sizes[mpi_rank],
		   MPI_DOUBLE,
		   predictions,
		   residual_sizes,
		   residual_offsets,
		   MPI_DOUBLE,
		   communicator);

    double likelihood = 0.0;
    
    if (mpi_rank == 0) {

      double hvalue = lambda_scale;
      
      if (tdtwave2d_compute_likelihood_(&nobservations,
					&hvalue,
					predictions,
					residual,
					unused,
					&likelihood,
					&log_normalization) < 0) {
	throw TDTWAVE2DEXCEPTION("Failed to compute likelihood");
      }

    }

    MPI_Bcast(&likelihood, 1, MPI_DOUBLE, 0, communicator);
    MPI_Bcast(&log_normalization, 1, MPI_DOUBLE, 0, communicator);

    return likelihood;
    
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
