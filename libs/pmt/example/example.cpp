#include <iostream>
#include <memory>

#include <pmt.h>

int main(int argc, char* argv[]) {
  std::unique_ptr<pmt::PMT> sensor = pmt::nvml::NVML::Create();
  auto start = sensor->Read();
  std::this_thread::sleep_for(
      std::chrono::milliseconds(sensor->GetMeasurementInterval()));
  auto end = sensor->Read();
  std::cout << "Runtime: " << pmt::PMT::seconds(start, end) << " s"
            << std::endl;
  std::cout << "Joules: " << pmt::PMT::joules(start, end) << " J" << std::endl;
  std::cout << "Watt: " << pmt::PMT::watts(start, end) << " W" << std::endl;

  return 0;
}