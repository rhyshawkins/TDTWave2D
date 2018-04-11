//
//    TDTWave2D : General framework for inversion using the
//    trans-dimensional tree approach with a wavelet parameterisation. See
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

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

		 
