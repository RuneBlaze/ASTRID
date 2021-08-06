#include "DistanceMethods.hpp"
#include "glog/logging.h"
#include <iostream>

std::string RapidNJ(TaxonSet &ts, DistanceMatrix &dm) {

  LOG(ERROR) << "RapidNJ is not supported on Windows" << std::endl;
  return "RapidNJ is not supported on Windows";
}