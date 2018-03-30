
#include <stdio.h>
#include <stdarg.h>

#include "tdt2dwaveexception.hpp"

extern "C" {
  #include "slog.h"
};

tdt2dwaveexception::tdt2dwaveexception(const char *srcfile,
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

tdt2dwaveexception::~tdt2dwaveexception()
{
}
