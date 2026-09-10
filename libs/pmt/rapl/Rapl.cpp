#include "Rapl.h"
#include "RaplImpl.h"

namespace pmt::rapl {

std::unique_ptr<Rapl> Rapl::Create() { return std::make_unique<RaplImpl>(); }

}  // end namespace pmt::rapl
