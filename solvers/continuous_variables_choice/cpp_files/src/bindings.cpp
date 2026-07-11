#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "cnrpa.hpp"
#include "gaco_cabgnrpa.hpp"

namespace py = pybind11;

// ---------------------------------------------------------------------------
// Shared evaluator: calls py_udp.fitness(x) from C++
// ---------------------------------------------------------------------------
static Evaluator make_evaluator(py::object py_udp) {
    return [py_udp](const std::vector<double>& x) -> double {
        py::gil_scoped_acquire acquire;
        py::object result = py_udp.attr("fitness")(x);
        return result.cast<std::vector<double>>()[0];
    };
}

// ---------------------------------------------------------------------------
// cnrpa  (unchanged)
// ---------------------------------------------------------------------------
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

    Evaluator evaluator = make_evaluator(py_evaluator);

    CNRPAResult result;
    {
        py::gil_scoped_release release;
        result = cnrpa(
            evaluator, level, n_policies, bounds,
            learning_rate, timeout,
            unfeasibility_value, gaussian_kernel_threshold, seed);
    }

    py::dict out;
    out["best_values_sequence"] = result.best_values_sequence;
    out["best_states_sequence"] = result.best_states_sequence;
    out["best_value"]           = result.best_value;
    out["best_values_list"]     = result.best_values_list;
    out["time_list"]            = result.time_list;
    return out;
}

// ---------------------------------------------------------------------------
// gaco_cabgnrpa  (now using pagmo directly — no pg.algorithm/pg.population)
// ---------------------------------------------------------------------------
py::dict py_gaco_cabgnrpa(
    py::object py_udp,             // pykep UDP object
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
    // Extract bounds as pagmo::vector_double
    pagmo::vector_double lb, ub;
    Bounds               bounds;

    for (auto item : py_bounds) {
        auto pair = item.cast<std::vector<double>>();
        bounds.emplace_back(pair[0], pair[1]);
        lb.push_back(pair[0]);
        ub.push_back(pair[1]);
    }

    Evaluator evaluator = make_evaluator(py_udp);

    GACOCNRPAResult result;
    {
        py::gil_scoped_release release;
        result = gaco_cabgnrpa(
            evaluator, py_udp, lb, ub,
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

// ---------------------------------------------------------------------------
// Module
// ---------------------------------------------------------------------------
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
        py::arg("seed")                      = 0u,
        "Run CNRPA in C++ with Python evaluator callback.");

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
        py::arg("random_seed")               = 42u,
        R"doc(
Run GACO-CABGNRPA entirely in C++ using pagmo natively.

GACO runs as pure C++ pagmo — no Python roundtrip for evolve().
Only fitness() callbacks re-enter Python (unavoidable: evaluator is pykep).

Parameters
----------
udp          : pykep UDP object (e.g. pk.trajopt.gym.cassini1)
level        : nesting depth
n_policies   : rollouts per level
zeta         : std_factor multiplier for fit_gaussian_from_density
bounds       : list of [lo, hi] pairs
tau          : mixing weight between policy and GACO bias
kernel_size  : GACO archive size (default 63)
n_generations: GACO generations per evolve() call (default 10)
elitism_factor: GACO q parameter (default 1.0)
random_seed  : seed for both GACO and CNRPA RNG
        )doc");

    py::class_<GaussianKernel>(m, "GaussianKernel")
        .def(py::init<const std::vector<double>&, double>())
        .def("pdf", [](const GaussianKernel& k, const std::vector<double>& x) {
            return k.pdf(x);
        });
}
