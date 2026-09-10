#ifndef LIKWID_PMT_H
#define LIKWID_PMT_H

#include <memory>
#include <string_view>

#include "common/PMT.h"

namespace pmt {
namespace likwid {
class Likwid : public PMT {
 public:
  constexpr static inline std::string_view name = "likwid";

  static std::unique_ptr<Likwid> Create(
      std::string event_group_name = default_event_group());

  static std::string default_event_group() { return "ENERGY"; }
};
}  // end namespace likwid
}  // end namespace pmt

#endif
