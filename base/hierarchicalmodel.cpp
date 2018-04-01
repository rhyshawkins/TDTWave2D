
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
