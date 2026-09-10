#ifndef PMT_ROCM_H_
#define PMT_ROCM_H_

#include <memory>
#include <string_view>

#include "common/PMT.h"

namespace pmt::rocm {
class ROCM : public PMT {
 public:
  constexpr static inline std::string_view name = "rocm";
  static std::unique_ptr<ROCM> Create(int device_number = 0);
};
}  // end namespace pmt::rocm

#endif  // PMT_ROCM_H_
