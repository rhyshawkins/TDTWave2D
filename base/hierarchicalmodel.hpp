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
#ifndef hierarchicalmodel_hpp
#define hierarchicalmodel_hpp

#include <ostream>
#include <istream>

class hierarchicalmodel {
public:

  hierarchicalmodel(int nhierarchical);
  virtual ~hierarchicalmodel();

  virtual int get_nhierarchical() const;
  virtual void set(int i, double v);
  virtual double get(int i) const;

  virtual bool write(std::ostream &s) const;
  virtual bool read(std::istream &s);

protected:

  int nhierarchical;
  double *hierarchical;
  
};

class singlescaling_hierarchicalmodel : public hierarchicalmodel {
public:

  singlescaling_hierarchicalmodel(double lambda = 1.0);
  ~singlescaling_hierarchicalmodel();

};

#endif // hierarchicalmodel_hpp
  
