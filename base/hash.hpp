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
#ifndef hash_hpp
#define hash_hpp

#include <openssl/md5.h>

//
// Requires -lssl -lcrypto
//

class hash {
public:

  hash(const double *v,
       size_t n);

  void compute(const double *v,
	       size_t n);

  bool operator==(const hash &rhs);

private:
  
  unsigned char c[MD5_DIGEST_LENGTH];
  
};

#endif // hash_hpp

