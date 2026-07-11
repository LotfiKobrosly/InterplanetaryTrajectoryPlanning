#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pagmo/types.hpp>

#include <vector>
#include <utility>
#include <memory>

namespace py = pybind11;

// ---------------------------------------------------------------------------
// PythonEvaluatorUDP
//
// pagmo User Defined Problem wrapping a Python pykep evaluator.
// pagmo requires UDP to be CopyConstructible — we satisfy this by storing
// the Python object via shared_ptr (shallow copy, both share the same object).
// This is safe since the Python object is read-only from GACO's perspective
// (GACO only calls fitness(), never modifies the evaluator).
// ---------------------------------------------------------------------------
class PythonEvaluatorUDP {
public:
    // Default constructor required by pagmo's is_udp type trait
    PythonEvaluatorUDP() = default;

    PythonEvaluatorUDP(
        py::object                   py_udp,
        const pagmo::vector_double&  lb,
        const pagmo::vector_double&  ub)
        : py_udp_(std::make_shared<py::object>(std::move(py_udp)))
        , lb_(lb)
        , ub_(ub)
    {}

    // Copy constructor: shared_ptr copy — both instances share py_udp_
    PythonEvaluatorUDP(const PythonEvaluatorUDP&)            = default;
    PythonEvaluatorUDP(PythonEvaluatorUDP&&)                 = default;
    PythonEvaluatorUDP& operator=(const PythonEvaluatorUDP&) = default;
    PythonEvaluatorUDP& operator=(PythonEvaluatorUDP&&)      = default;

    // pagmo UDP interface: fitness
    pagmo::vector_double fitness(const pagmo::vector_double& x) const {
        if (!py_udp_)
            throw std::runtime_error(
                "PythonEvaluatorUDP: fitness() called on default-constructed instance");
        py::gil_scoped_acquire acquire;
        py::object result = py_udp_->attr("fitness")(x);
        return result.cast<pagmo::vector_double>();
    }

    // pagmo UDP interface: bounds
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const {
        return {lb_, ub_};
    }

private:
    std::shared_ptr<py::object> py_udp_;
    pagmo::vector_double        lb_;
    pagmo::vector_double        ub_;
};
