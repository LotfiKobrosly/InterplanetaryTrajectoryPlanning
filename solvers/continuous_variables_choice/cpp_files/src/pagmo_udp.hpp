#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pagmo/types.hpp>

#include <vector>
#include <utility>
#include <memory>
#include <cstdio>

namespace py = pybind11;

class PythonEvaluatorUDP {
public:
    PythonEvaluatorUDP() = default;

    PythonEvaluatorUDP(
        py::object                   py_udp,
        const pagmo::vector_double&  lb,
        const pagmo::vector_double&  ub)
        : py_udp_(std::make_shared<py::object>(std::move(py_udp)))
        , lb_(lb)
        , ub_(ub)
    {}

    PythonEvaluatorUDP(const PythonEvaluatorUDP&)            = default;
    PythonEvaluatorUDP(PythonEvaluatorUDP&&)                 = default;
    PythonEvaluatorUDP& operator=(const PythonEvaluatorUDP&) = default;
    PythonEvaluatorUDP& operator=(PythonEvaluatorUDP&&)      = default;

    pagmo::vector_double fitness(const pagmo::vector_double& x) const {
        if (!py_udp_) {
            fprintf(stderr, "[UDP] WARNING: default-constructed fitness called\n");
            fflush(stderr);
            return {1e10};
        }
        py::gil_scoped_acquire acquire;

        static int call_count = 0;
        ++call_count;
        if (call_count <= 3) {
            fprintf(stderr, "[UDP] gaco fitness call #%d x[0]=%.4f\n",
                    call_count, x.empty() ? 0.0 : x[0]);
            fflush(stderr);
        }

        py::object result = py_udp_->attr("fitness")(x);
        auto r = result.cast<pagmo::vector_double>();

        if (call_count <= 3) {
            fprintf(stderr, "[UDP] gaco fitness result=%.4f\n",
                    r.empty() ? 0.0 : r[0]);
            fflush(stderr);
        }
        return r;
    }

    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const {
        return {lb_, ub_};
    }

private:
    std::shared_ptr<py::object> py_udp_;
    pagmo::vector_double        lb_;
    pagmo::vector_double        ub_;
};