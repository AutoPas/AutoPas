#ifndef PMT_RAPL_H_
#define PMT_RAPL_H_

#include <memory>
#include <string_view>

#include "common/PMT.h"

namespace pmt::rapl {
class Rapl : public PMT {
 public:
  constexpr static inline std::string_view name = "rapl";
  static std::unique_ptr<Rapl> Create();
};
}  // end namespace pmt::rapl

#endif  // PMT_RAPL_H_
