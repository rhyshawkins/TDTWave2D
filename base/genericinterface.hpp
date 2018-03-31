#pragma once
#ifndef genericinterface_hpp
#define genericinterface_hpp

extern "C" {

  
  //
  // Load data and return the number of observations
  // on success and -1 on error
  //
  int tdtwave2d_loaddata_(int *n,
			  const char *filename,
			  int *width,
			  int *height);

  //
  // For a single observation, compute the prediction given
  // model image
  //
  int tdtwave2d_compute_prediction_(int *observation,
				    int *width,
				    int *height,
				    const double *image,
				    double *unused,
				    double *prediction);
  //
  // For the observations, given predictions, compute residuals, likelihood
  // and norm
  //
  int tdtwave2d_compute_likelihood_(int *nobservation,
				    double *hierarchical,
				    double *predictions,
				    double *residuals,
				    double *unused,
				    double *like,
				    double *norm);

  //
  // Used for making synthetic datasets, save data in correct format
  // using predictions to overwrite observations
  //
  int tdtwave2d_savedata_(int *n,
		     const char *filename,
		     double *hierarchical,
		     int *nobservations,
		     double *predictions);

  //
  // Finish (save state for restarting)
  //
  int tdtwave2d_finish_(const char *path);

  //
  // Starting (pre load)
  //
  int tdtwave2d_start_(const char *path);

};

#endif // genericinterface_hpp

		 
