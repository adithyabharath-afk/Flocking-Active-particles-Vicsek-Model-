#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "simulation_flock.h"

namespace py = pybind11;

PYBIND11_MODULE(flocksim, m) {
    py::class_<simulation_flock>(m, "Simulation_flock")
        .def(py::init<double, double, int>(),
             py::arg("length"), py::arg("cutoff"), py::arg("num_particles"))
        .def("makegrid", &simulation_flock::makegrid, "Build the neighbor list grid",
             py::arg("coordinates"))
        .def("force", &simulation_flock::force, "Compute 2D forces using the grid",
             py::arg("coordinates"), py::arg("orientations"));
}