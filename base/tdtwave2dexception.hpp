#pragma once
#ifndef tdtwave2dexception_hpp
#define tdtwave2dexception_hpp

#include <exception>

#define TDTWAVE2DEXCEPTION(fmt, ...) tdtwave2dexception(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

class tdtwave2dexception : public std::exception {
public:

  
  tdtwave2dexception(const char *srcfile,
		     const char *function,
		     int lineno,
		     const char *fmt, ...);
  ~tdtwave2dexception();
  
};

#endif // tdtwave2dexception_hpp
