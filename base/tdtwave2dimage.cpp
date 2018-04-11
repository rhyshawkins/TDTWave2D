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

