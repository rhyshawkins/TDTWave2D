#include <vector>
#include <stdio.h>
#include <math.h>

#include "genericinterface.hpp"

#include "linearweights.hpp"

struct observation {
  double value;
  double sigma;

  std::vector<double> lon;
  std::vector<double> lat;
  LinearWeights weights;
};

struct model_bounds {
  double minlon;
  double minlat;
  double maxlat;
  double maxlon;
};

static constexpr double RADIUS = 6371.0;
  
static struct model_bounds bounds;
static std::vector<struct observation> obs;

static double mean_offset = 0.0;

extern "C" {

  int tdtwave2d_loaddata_(int *n,
			  const char *filename,
			  int *width,
			  int *height)
  {
    printf("here\n");
    
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      return -1;
    }

    if (fscanf(fp, "%lf %lf %lf %lf\n",
	       &bounds.minlon,
	       &bounds.maxlon,
	       &bounds.minlat,
	       &bounds.maxlat) != 4) {
      return -1;
    }
	
    int nobs;
    if (fscanf(fp, "%d\n", &nobs) != 1) {
      return -1;
    }

    obs.resize(nobs);
    mean_offset = 0.0;
    for (int i = 0; i < nobs; i ++) {

      int npoints;
      if (fscanf(fp, "%lf %lf %d\n",
		 &obs[i].value,
		 &obs[i].sigma,
		 &npoints) != 3) {
	return -1;
      }

      if (npoints <= 1) {
	fprintf(stderr, "error: need at least two points for a path\n");
	return -1;
      }

      obs[i].lon.resize(npoints);
      obs[i].lat.resize(npoints);
      
      for (int j = 0; j < npoints; j ++) {
	if (fscanf(fp, "%lf %lf\n",
		   &obs[i].lon[j],
		   &obs[i].lat[j]) != 2) {
	  return -1;
	}
      }

      obs[i].weights.clear();
      for (int j = 1; j < npoints; j ++) {
	if (!obs[i].weights.compute_weights(bounds.minlon, bounds.maxlon,
					    bounds.minlat, bounds.maxlat,
					    *width, *height,
					    RADIUS,
					    obs[i].lon[j - 1], obs[i].lat[j - 1],
					    obs[i].lon[j], obs[i].lat[j])) {
	  fprintf(stderr, "error: failed to compute weights for path\n");
	  return -1;
	}
      }

      if (obs[i].value > 0.0) {
	mean_offset += obs[i].weights.total_weight()/obs[i].value;
      }
    }

    if (mean_offset > 0.0) {
      mean_offset /= (double)nobs;
    }

    printf("Mean velocity: %10.6f\n", mean_offset);
    
    fclose(fp);
    return nobs;
  }

  int tdtwave2d_compute_prediction_(int *observation,
				    int *width,
				    int *height,
				    const double *image,
				    double *unused,
				    double *prediction)
  {
    if ((*observation) < (int)obs.size()) {
      
      prediction[0] = obs[*observation].weights.evaluate_velocityfield(mean_offset, image);
			    
      return 0;
    }

    return -1;
  }

  int tdtwave2d_compute_likelihood_(int *nobservation,
				    double *hierarchical,
				    double *predictions,
				    double *residuals,
				    double *unused,
				    double *_like,
				    double *_norm)

  {
    double sum = 0.0;
    double norm = -0.5*log(2.0*M_PI);
    for (int i = 0; i < (*nobservation); i ++) {
      double res = predictions[i] - obs[i].value;
      residuals[i] = res;
      double n = hierarchical[0] * obs[i].sigma;

      sum += (res*res)/(2.0*n*n);
      norm += log(n);
    }

    *_like = sum;
    *_norm = norm;

    return 0;
  }

  int tdtwave2d_savedata_(int *n,
			  const char *filename,
			  double *noiselevel,
			  int *nobservations,
			  double *predictions)
  {
    FILE *fp;

    fp = fopen(filename, "w");
    if (fp == NULL) {
      return -1;
    }

    fprintf(fp, "%15.9f %15.9f %15.9f %15.9f\n",
	    bounds.minlon,
	    bounds.maxlon,
	    bounds.minlat,
	    bounds.maxlat);
	
    fprintf(fp, "%d\n", (int)obs.size());
    int i = 0;
    for (auto &o : obs) {
      
      fprintf(fp, "%16.9f %16.9f %d\n",
	      predictions[i], noiselevel[0], (int)o.lon.size());

      for (int j = 0; j < (int)o.lon.size(); j ++) {
	fprintf(fp, "%16.9f %16.9f\n", o.lon[j], o.lat[j]);
      }

      i ++;
    }

    fclose(fp);
    return 0;
  }

  
}

