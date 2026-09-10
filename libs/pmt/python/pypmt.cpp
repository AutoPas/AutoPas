#include <pybind11/pybind11.h>

#include <pmt.h>

namespace py = pybind11;

PYBIND11_MODULE(pmt, m) {
  m.doc() = "libpmt python bindings";

  m.def(
      "create", [](py::str name) { return pmt::Create(name, ""); },
      "Create PMT instance");

  m.def(
      "create",
      [](py::str name, py::object argument) {
        return pmt::Create(name, py::str(argument));
      },
      "Create PMT instance");

  py::class_<pmt::PMT>(m, "PMT")
      .def("seconds",
           py::overload_cast<const pmt::State&, const pmt::State&>(
               &pmt::PMT::seconds),
           "Get elapsed time")
      .def("joules", &pmt::PMT::joules, "Get energy consumption")
      .def("watts", &pmt::PMT::watts, "Get average power consumption")
      .def("read", &pmt::PMT::Read)
      .def("startDump", &pmt::PMT::StartDump)
      .def("stopDump", &pmt::PMT::StopDump);

  py::class_<pmt::State>(m, "State")
      .def("timestamp", &pmt::State::timestamp, "Get timestamp")
      .def("watts", &pmt::State::watts, py::arg("index"),
           "Get instantenous power consumption for the specified measurement")
      .def("name", &pmt::State::name, py::arg("index"),
           "Get name for the specified measurement")
      .def("nr_measurements", &pmt::State::NrMeasurements,
           "Get number of distinct measurements");
}
