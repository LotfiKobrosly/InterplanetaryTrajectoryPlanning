#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "cnrpa.hpp"
#include "gaco_cabgnrpa.hpp"
#include "pagmo_udp.hpp"

namespace py = pybind11;

// ---------------------------------------------------------------------------
// make_evaluator
// Passes x as a numpy array to avoid std::vector->list conversion cost.
// This reduces pybind11 overhead from ~44us to ~5us per call.
// GIL is held throughout CNRPA — no acquire/release needed per call.
// ---------------------------------------------------------------------------
static Evaluator make_evaluator_gil_held(py::object py_udp) {
    return [py_udp](const std::vector<double>& x) -> double {
        // Wrap x as numpy array without copying — zero allocation
        py::array_t<double> arr(
            {(py::ssize_t)x.size()},   // shape
            {sizeof(double)},           // strides
            x.data(),                   // pointer to data (no copy)
            py::none()                  // base object (none = C++ owns data)
        );
        py::object result = py_udp.attr("fitness")(arr);
        return result.cast<std::vector<double>>()[0];
    };
}

// For gaco_cabgnrpa: GIL released for pagmo, must reacquire per call
static Evaluator make_evaluator_gil_released(py::object py_udp) {
    return [py_udp](const std::vector<double>& x) -> double {
        py::gil_scoped_acquire acquire;
        py::array_t<double> arr(
            {(py::ssize_t)x.size()},
            {sizeof(double)},
            x.data(),
            py::none()
        );
        py::object result = py_udp.attr("fitness")(arr);
        return result.cast<std::vector<double>>()[0];
    };
}

py::dict py_cnrpa(
    py::object py_evaluator,
    int        level,
    int        n_policies,
    py::list   py_bounds,
    double     learning_rate,
    double     timeout,
    double     unfeasibility_value,
    double     gaussian_kernel_threshold,
    uint64_t   seed)
{
    Bounds bounds;
    for (auto item : py_bounds) {
        auto pair = item.cast<std::vector<double>>();
        bounds.emplace_back(pair[0], pair[1]);
    }

    Evaluator evaluator = make_evaluator_gil_held(py_evaluator);

    // Do NOT release GIL — fitness calls back into Python anyway
    CNRPAResult result = cnrpa(
        evaluator, level, n_policies, bounds,
        learning_rate, timeout,
        unfeasibility_value, gaussian_kernel_threshold, seed);

    py::dict out;
    out["best_values_sequence"] = result.best_values_sequence;
    out["best_states_sequence"] = result.best_states_sequence;
    out["best_value"]           = result.best_value;
    out["best_values_list"]     = result.best_values_list;
    out["time_list"]            = result.time_list;
    return out;
}

py::dict py_gaco_cabgnrpa(
    py::object py_udp,
    int        level,
    int        n_policies,
    double     zeta,
    py::list   py_bounds,
    double     learning_rate,
    double     timeout,
    double     tau,
    double     unfeasibility_value,
    double     gaussian_kernel_threshold,
    unsigned   kernel_size,
    unsigned   n_generations,
    double     elitism_factor,
    unsigned   random_seed)
{
    pagmo::vector_double lb, ub;
    Bounds               bounds;
    for (auto item : py_bounds) {
        auto pair = item.cast<std::vector<double>>();
        bounds.emplace_back(pair[0], pair[1]);
        lb.push_back(pair[0]);
        ub.push_back(pair[1]);
    }

    // Build pagmo::problem WITH GIL held
    pagmo::problem prob{PythonEvaluatorUDP(py_udp, lb, ub)};
    Evaluator evaluator = make_evaluator_gil_released(py_udp);

    GACOCNRPAResult result;
    {
        py::gil_scoped_release release;
        result = gaco_cabgnrpa(
            evaluator, prob, lb, ub,
            level, n_policies, zeta, bounds,
            learning_rate, timeout, tau,
            unfeasibility_value, gaussian_kernel_threshold,
            kernel_size, n_generations, elitism_factor, random_seed);
    }

    py::dict out;
    out["best_values_sequence"] = result.best_values_sequence;
    out["best_states_sequence"] = result.best_states_sequence;
    out["best_value"]           = result.best_value;
    out["best_values_list"]     = result.best_values_list;
    out["time_list"]            = result.time_list;
    return out;
}

PYBIND11_MODULE(cnrpa_cpp, m) {
    m.doc() = "C++ CNRPA and GACO-CABGNRPA using pagmo natively";

    m.def("cnrpa", &py_cnrpa,
        py::arg("evaluator"),
        py::arg("level"),
        py::arg("n_policies"),
        py::arg("bounds"),
        py::arg("learning_rate")             = 0.01,
        py::arg("timeout")                   = 10.0,
        py::arg("unfeasibility_value")       = 1e+10,
        py::arg("gaussian_kernel_threshold") = 1e-6,
        py::arg("seed")                      = 0u);

    m.def("gaco_cabgnrpa", &py_gaco_cabgnrpa,
        py::arg("udp"),
        py::arg("level"),
        py::arg("n_policies"),
        py::arg("zeta"),
        py::arg("bounds"),
        py::arg("learning_rate")             = 0.01,
        py::arg("timeout")                   = 10.0,
        py::arg("tau")                       = 10.0,
        py::arg("unfeasibility_value")       = 1e+10,
        py::arg("gaussian_kernel_threshold") = 1e-6,
        py::arg("kernel_size")               = 63u,
        py::arg("n_generations")             = 10u,
        py::arg("elitism_factor")            = 1.0,
        py::arg("random_seed")               = 42u);

    py::class_<GaussianKernel>(m, "GaussianKernel")
        .def(py::init<const std::vector<double>&, double>())
        .def("pdf", [](const GaussianKernel& k, const std::vector<double>& x) {
            return k.pdf(x);
        });
}