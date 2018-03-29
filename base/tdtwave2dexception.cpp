
#include <stdio.h>
#include <stdarg.h>

#include "aemexception.hpp"

extern "C" {
  #include "slog.h"
};

aemexception::aemexception(const char *srcfile,
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

aemexception::~aemexception()
{
}
