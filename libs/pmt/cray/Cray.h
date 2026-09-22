#ifndef CRAY_PMT_H
#define CRAY_PMT_H

#include <cstddef>
#include <memory>
#include <string_view>

#include "common/PMT.h"

namespace pmt {
namespace cray {

struct CrayMeasurement {
  std::string name;
  size_t watt;
};

class Cray : public PMT {
 public:
  constexpr static inline std::string_view name = "cray";
  static std::unique_ptr<Cray> Create();
};
}  // end namespace cray
}  // end namespace pmt

#endif
