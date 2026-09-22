#ifndef PMT_NVIDIA_H_
#define PMT_NVIDIA_H_

#include <memory>
#include <string_view>

#include "common/PMT.h"

namespace pmt::nvidia {
class NVIDIA : public PMT {
 public:
  constexpr static inline std::string_view name = "nvidia";
  static std::unique_ptr<PMT> Create(int device_number = 0);
};
}  // end namespace pmt::nvidia

#endif  // PMT_NVIDIA_H_
