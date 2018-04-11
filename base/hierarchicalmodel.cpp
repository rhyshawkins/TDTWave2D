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

#include "hierarchicalmodel.hpp"

#include "tdtwave2dexception.hpp"

hierarchicalmodel::hierarchicalmodel(int _nhierarchical) :
  nhierarchical(_nhierarchical),
  hierarchical(new double[_nhierarchical])
{
}

hierarchicalmodel::~hierarchicalmodel()
{
  delete [] hierarchical;
}

int
hierarchicalmodel::get_nhierarchical() const
{
  return nhierarchical;
}

void
hierarchicalmodel::set(int i, double v)
{
  if (i < 0 || i >= nhierarchical) {
    throw TDTWAVE2DEXCEPTION("Index out of range %d (%d)\n", i, nhierarchical);
  }

  hierarchical[i] = v;
}

double
hierarchicalmodel::get(int i) const
{
  if (i < 0 || i >= nhierarchical) {
    throw TDTWAVE2DEXCEPTION("Index out of range %d (%d)\n", i, nhierarchical);
  }

  return hierarchical[i];
}

bool
hierarchicalmodel::write(std::ostream &s) const
{
  s.write((char*)&nhierarchical, sizeof(int));
  
  for (int i = 0; i < nhierarchical; i ++) {
    s.write((char*)(&hierarchical[i]), sizeof(double));
  }

  return true;
}

bool
hierarchicalmodel::read(std::istream &s)
{
  int n;

  s.read((char*)&n, sizeof(int));
  if (n != nhierarchical) {
    return false;
  }

  for (int i = 0; i < nhierarchical; i ++) {
    s.read((char*)(&hierarchical[i]), sizeof(double));
  }

  return true;
}

singlescaling_hierarchicalmodel::singlescaling_hierarchicalmodel(double lambda) :
  hierarchicalmodel(1)
{
  set(0, lambda);
}

singlescaling_hierarchicalmodel::~singlescaling_hierarchicalmodel()
{
}
