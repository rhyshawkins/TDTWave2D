#pragma once
#ifndef tdt2dwaveexception_hpp
#define tdt2dwaveexception_hpp

#include <exception>

#define TDT2DWAVEEXCEPTION(fmt, ...) tdt2dwaveexception(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

class tdt2dwaveexception : public std::exception {
public:

  
  tdt2dwaveexception(const char *srcfile,
		     const char *function,
		     int lineno,
		     const char *fmt, ...);
  ~tdt2dwaveexception();
  
};

#endif // tdt2dwaveexception_hpp
