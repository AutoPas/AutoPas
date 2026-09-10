#include <utility>

#include "Likwid.h"
#include "LikwidImpl.h"

namespace pmt::likwid {

std::unique_ptr<Likwid> Likwid::Create(std::string event_group_name) {
  return std::make_unique<LikwidImpl>(std::move(event_group_name));
}

}  // namespace pmt::likwid
