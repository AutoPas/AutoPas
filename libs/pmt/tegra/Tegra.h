#ifndef TegraImplPMT_H
#define TegraImplPMT_H

#include <memory>
#include <string_view>

#include "common/PMT.h"

namespace pmt::tegra {
class Tegra : public PMT {
 public:
  constexpr static inline std::string_view name = "tegra";
  static std::unique_ptr<Tegra> Create();
};
}  // end namespace pmt::tegra

#endif
