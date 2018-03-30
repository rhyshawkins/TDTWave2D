
#include <stdio.h>
#include <stdarg.h>

#include "tdtwave2dexception.hpp"

extern "C" {
  #include "slog.h"
};

tdtwave2dexception::tdtwave2dexception(const char *srcfile,
			   const char *function,
			   int lineno,
			   const char *fmt, ...)
{
  va_list ap;
  
  va_start(ap, fmt);
  vslog(SLOG_ERROR,
	srcfile,
	function,
	lineno,
	fmt,
	ap);
  va_end(ap);
  
}

tdtwave2dexception::~tdtwave2dexception()
{
}
