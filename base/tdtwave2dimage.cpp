
#include <stdio.h>
#include <math.h>

#include "tdtwave2dimage.hpp"
#include "tdtwave2dexception.hpp"



tdtwave2dimage::tdtwave2dimage() :
  rows(-1),
  columns(-1),
  image(nullptr)
{
}

tdtwave2dimage::tdtwave2dimage(int _rows, int _columns, double const_image) :
  rows(_rows),
  columns(_columns),
  image(new double[rows * columns])
{
  int s = rows * columns;
  for (int i = 0; i < s; i ++) {
    image[i] = const_image;
  }
}

tdtwave2dimage::~tdtwave2dimage()
{
  delete [] image;
}

bool
tdtwave2dimage::load(const char *filename)
{
  int nrows, ncols;
  
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

  image = new double[rows * columns];

  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {

      if (fscanf(fp, "%lf", &image[j * columns + i]) != 1) {
	return false;
      }
    }
  }
  
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

