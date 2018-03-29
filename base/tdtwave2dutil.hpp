#pragma once
#ifndef aemutil_hpp

#include <string>
#include <vector>

std::string mkfilename(const char *prefix, const char *file);

std::string mkfilenamerank(const char *prefix, const char *file, int rank);

std::string mkformatstring(const char *fmt, ...);

std::string stripwhitespaceandquotes(const char *s);

bool loadhierarchicallambda(const char *filename, std::vector<double> &lambda);


#endif // aemutil.hpp
