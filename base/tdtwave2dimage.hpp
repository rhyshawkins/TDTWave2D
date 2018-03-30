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
  
