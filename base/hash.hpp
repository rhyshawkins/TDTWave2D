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

