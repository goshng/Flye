//(c) 2013-2016 by Authors
// This file is a part of Ragout program.
// Released under the BSD license (see LICENSE file)

#pragma once

#include <ctime>
#include <iostream>
#include <iterator>
#include <sstream>
#include <utility>
#include <vector>

inline std::vector<std::string> splitString(const std::string &s, char delim) {
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim))
    elems.push_back(item);
  return elems;
}
