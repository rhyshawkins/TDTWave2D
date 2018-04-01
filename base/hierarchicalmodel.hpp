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
  
