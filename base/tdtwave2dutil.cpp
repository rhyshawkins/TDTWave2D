
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "aemutil.hpp"

std::string mkfilename(const char *prefix, const char *file)
{
  if (prefix == nullptr) {
    return std::string(file);
  } else {
    return std::string(prefix) + file;
  }
}

std::string mkfilenamerank(const char *prefix, const char *file, int rank)
{
  char buffer[16];

  sprintf(buffer, "-%03d", rank);
  
  if (prefix == nullptr) {
    return std::string(file) + buffer;
  } else {
    return std::string(prefix) + file + buffer;
  }
}

std::string mkformatstring(const char *fmt, ...)
{
  static char *buffer = nullptr;
  static int buffer_size = -1;

  if (buffer == nullptr) {
    buffer_size = 512;
    buffer = new char[buffer_size];
  }

  va_list ap;
  int size;
  
  va_start(ap, fmt);
  size = vsnprintf(buffer, buffer_size, fmt, ap);
  while (size >= buffer_size) {
    delete [] buffer;
    buffer_size *= 2;
    buffer = new char[buffer_size];
    size = vsnprintf(buffer, buffer_size, fmt, ap);
  }
  va_end(ap);

  return std::string(buffer);
}

std::string stripwhitespaceandquotes(const char *s)
{
  int i = 0;
  int j = strlen(s);

  while (i < j && (isspace(s[i]) || s[i] == '"' || s[j] == '\'')) {
    i ++;
  }

  while (j > i && (isspace(s[j]) || s[j] == '"' || s[j] == '\'')) {
    j --;
  }

  return std::string(s + i, j - i - 1);
}

bool loadhierarchicallambda(const char *filename, std::vector<double> &lambda)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    return false;
  }

  while (!feof(fp)) {
    double d;
    if (fscanf(fp, "%lf\n", &d) != 1) {
      if (!feof(fp)) {
	return false;
      } else {
	break;
      }
    } else {
      lambda.push_back(d);
    }
  }

  fclose(fp);
  return true;
}

