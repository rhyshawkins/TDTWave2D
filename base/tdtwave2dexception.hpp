#pragma once
#ifndef aemexception_hpp
#define aemexception_hpp

#include <exception>

#define AEMEXCEPTION(fmt, ...) aemexception(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

class aemexception : public std::exception {
public:

  
  aemexception(const char *srcfile,
		const char *function,
		int lineno,
		const char *fmt, ...);
  ~aemexception();
  
};

#endif // aemexception_hpp
