#include "Tegra.h"
#include "TegraImpl.h"

namespace pmt::tegra {

std::unique_ptr<Tegra> Tegra::Create() { return std::make_unique<TegraImpl>(); }

}  // end namespace pmt::tegra
