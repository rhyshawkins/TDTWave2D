
#include "hash.hpp"

hash::hash(const double *v,
	   size_t n)
{
  compute(v, n);
}

void
hash::compute(const double *v,
	      size_t n)
{
  MD5_CTX mdContext;

  MD5_Init(&mdContext);

  MD5_Update(&mdContext, v, n * sizeof(double));

  MD5_Final(c, &mdContext);
}

bool
hash::operator==(const hash &rhs)
{
  for (int i = 0; i < MD5_DIGEST_LENGTH; i ++) {

    if (c[i] != rhs.c[i]) {
      return false;
    }

  }

  return true;
}

