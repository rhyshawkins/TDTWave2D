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
#ifndef tdtwave2dimage_hpp
#define tdtwave2dimage_hpp

#include <vector>

class tdtwave2dimage {
public:

  tdtwave2dimage();
  tdtwave2dimage(int rows, int columns, double const_image = 0.0);
  ~tdtwave2dimage();

  bool load(const char *filename);
  bool save(const char *filename);

  bool save_image(const char *filename);

  int rows;
  int columns;
  
  double *image;

};

#endif // tdtwave2dimage_hpp
  
