
#include <stdio.h>
#include <math.h>

#include "tdtwave2dimage.hpp"
#include "tdtwave2dexception.hpp"



tdtwave2dimage::tdtwave2dimage() :
  rows(-1),
  columns(-1),
  depth(0.0),
  image(nullptr)
{
}

tdtwave2dimage::tdtwave2dimage(int _rows, int _columns, double _depth, double const_image) :
  rows(_rows),
  columns(_columns),
  depth(_depth),
  image(new double[rows * columns])
{
  int s = rows * columns;
  for (int i = 0; i < s; i ++) {
    image[i] = const_image;
  }

  update_layer_thickness();
}

tdtwave2dimage::~tdtwave2dimage()
{
  delete [] image;
}

bool
tdtwave2dimage::load(const char *filename)
{
  int nrows, ncols;
  double ndepth;
  
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    return false;
  }

  if (fscanf(fp, "%d %d\n", &nrows, &ncols) != 2) {
    return false;
  }
  
  if (image != nullptr) {
    delete [] image;
  }

  rows = nrows;
  columns = ncols;
  depth = ndepth;

  image = new double[rows * columns];

  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {

      if (fscanf(fp, "%lf", &image[j * columns + i]) != 1) {
	return false;
      }
    }
  }

  update_layer_thickness();
  
  fclose(fp);
  return true;
}

bool
tdtwave2dimage::save(const char *filename)
{
  if (image == nullptr) {
    return false;
  }

  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    return false;
  }

  fprintf(fp, "%d %d\n", rows, columns);

  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {

      fprintf(fp, "%15.9f ", image[j * columns + i]);

    }

    fprintf(fp, "\n");
  }

  fclose(fp);

  return true;
}

bool
tdtwave2dimage::save_image(const char *filename)
{
  if (image == nullptr) {
    return false;
  }

  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    return false;
  }

  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {

      fprintf(fp, "%15.9f ", image[j * columns + i]);

    }

    fprintf(fp, "\n");
  }

  fclose(fp);

  return true;
}

