#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "cnrpa.hpp"

namespace py = pybind11;

// ---------------------------------------------------------------------------
// Helper: convert Python evaluator (any callable returning float) to
// C++ std::function<double(std::vector<double>)>
// ---------------------------------------------------------------------------
static Evaluator make_evaluator(py::object py_evaluator) {
    return [py_evaluator](const std::vector<double>& x) -> double {
        py::gil_scoped_acquire acquire;  // re-acquire GIL for Python call
        py::object result = py_evaluator.attr("fitness")(x);
        // fitness() returns a list/array with one element
        return result.cast<std::vector<double>>()[0];
    };
}

// ---------------------------------------------------------------------------
// Python-facing cnrpa() wrapper
// ---------------------------------------------------------------------------
py::dict py_cnrpa(
    py::object   py_evaluator,
    int          level,
    int          n_policies,
    py::list     py_bounds,       // list of [lo, hi] pairs
    double       learning_rate,
    double       timeout,
    double       unfeasibility_value,
    double       gaussian_kernel_threshold,
    uint64_t     seed)
{
    // Convert Python bounds [[lo, hi], ...] → C++ Bounds
    Bounds bounds;
    for (auto item : py_bounds) {
        auto pair = item.cast<std::vector<double>>();
        if (pair.size() != 2)
            throw std::invalid_argument("Each bound must be [lo, hi]");
        bounds.emplace_back(pair[0], pair[1]);
    }

    Evaluator evaluator = make_evaluator(py_evaluator);

    // Release GIL while running pure C++ code
    CNRPAResult result;
    {
        py::gil_scoped_release release;
        result = cnrpa(
            evaluator,
            level,
            n_policies,
            bounds,
            learning_rate,
            timeout,
            unfeasibility_value,
            gaussian_kernel_threshold,
            seed);
    }

    // Pack results into a Python dict
    py::dict out;
    out["best_values_sequence"] = result.best_values_sequence;
    out["best_states_sequence"] = result.best_states_sequence;
    out["best_value"]           = result.best_value;
    out["best_values_list"]     = result.best_values_list;
    out["time_list"]            = result.time_list;
    return out;
}

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------
PYBIND11_MODULE(cnrpa_cpp, m) {
    m.doc() = "C++ implementation of CNRPA for interplanetary trajectory optimization";

    m.def("cnrpa",
        &py_cnrpa,
        py::arg("evaluator"),
        py::arg("level"),
        py::arg("n_policies"),
        py::arg("bounds"),
        py::arg("learning_rate")             = 0.01,
        py::arg("timeout")                   = 10.0,
        py::arg("unfeasibility_value")       = 1e+10,
        py::arg("gaussian_kernel_threshold") = 1e-6,
        py::arg("seed")                      = 0,
        R"doc(
Run CNRPA (Continuous Nested Rollout Policy Adaptation) in C++,
calling back into a Python pykep evaluator for fitness evaluation.

Parameters
----------
evaluator               : object with a .fitness(x) method (e.g. pk.trajopt.mga instance)
level                   : nesting depth (0 = single playout)
n_policies              : number of rollouts per level
bounds                  : list of [lo, hi] pairs, one per decision variable
learning_rate           : policy adaptation step size (default 0.01)
timeout                 : wall-clock time limit in seconds (default 10.0)
unfeasibility_value     : fitness value that signals an infeasible solution (default 1e10)
gaussian_kernel_threshold: minimum kernel weight to include in policy (default 1e-6)
seed                    : RNG seed (0 = random device)

Returns
-------
dict with keys:
    best_values_sequence : list[float]  — best decision vector found
    best_states_sequence : list[float]  — corresponding normalized states
    best_value           : float        — best fitness value
    best_values_list     : list[float]  — history of best values over time
    time_list            : list[float]  — timestamps of improvements (seconds)
        )doc"
    );

    // Also expose GaussianKernel for testing/inspection from Python
    py::class_<GaussianKernel>(m, "GaussianKernel")
        .def(py::init<const std::vector<double>&, double>(),
             py::arg("center"), py::arg("sigma"))
        .def("pdf",
             [](const GaussianKernel& k, const std::vector<double>& x) {
                 return k.pdf(x);
             },
             py::arg("x"))
        .def("pdf_scalar",
             [](const GaussianKernel& k, double x) { return k.pdf(x); },
             py::arg("x"));
}
