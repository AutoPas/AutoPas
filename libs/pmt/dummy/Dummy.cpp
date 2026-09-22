#include <ext/alloc_traits.h>
#include <vector>

#include "Dummy.h"

namespace pmt {

class DummyImpl : public Dummy {
 private:
  virtual State GetState() override;

  virtual const char *GetDumpFilename() override { return nullptr; }
};

std::unique_ptr<Dummy> Dummy::Create() { return std::make_unique<DummyImpl>(); }

State DummyImpl::GetState() {
  State state;
  state.timestamp_ = PMT::GetTime();
  state.name_[0] = "none";
  state.joules_[0] = 0;
  state.watt_[0] = 0;
  return state;
}

}  // end namespace pmt
