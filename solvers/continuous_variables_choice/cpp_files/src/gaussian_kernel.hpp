#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

// ---------------------------------------------------------------------------
// GaussianKernel
// Mirrors the Python GaussianKernel class exactly:
//
//   norm_const = 1 / ((2 * pi * sigma) ^ (dim / 2))
//   pdf(x)     = norm_const * exp(-0.5 * ||x - mu||^2 / sigma^2)
//
// Note: this is an isotropic Gaussian (same sigma for all dimensions).
// The normalization constant includes sigma (not sigma^2) in the base,
// matching the Python implementation exactly.
// ---------------------------------------------------------------------------
class GaussianKernel {
public:
    GaussianKernel(const std::vector<double>& center, double sigma)
        : mu_(center), sigma_(sigma), dim_(center.size())
    {
        if (sigma <= 0.0)
            throw std::invalid_argument("GaussianKernel: sigma must be > 0");
        if (center.empty())
            throw std::invalid_argument("GaussianKernel: center must be non-empty");

        // Precompute normalization constant:
        // norm_const = 1 / ((2 * pi * sigma) ^ (dim / 2))
        double base   = 2.0 * M_PI * sigma_;
        norm_const_   = 1.0 / std::pow(base, static_cast<double>(dim_) / 2.0);
    }

    // Scalar overload: wraps x in a 1-D vector, matching Python's
    //   if isinstance(x, (float, int)): x = [x]
    GaussianKernel(double center_scalar, double sigma)
        : GaussianKernel(std::vector<double>{center_scalar}, sigma) {}

    // Evaluate pdf at point x (vector).
    // If x has fewer dimensions than the center, only the first x.size()
    // dimensions are used — matching numpy's silent truncation behavior
    // when np.dot() receives mismatched vectors (as in the original Python
    // where kernel center = states_sequence (N-D) but pdf receives a scalar).
    double pdf(const std::vector<double>& x) const {
        std::size_t effective_dim = std::min(dim_, x.size());

        if (effective_dim == 0)
            throw std::invalid_argument(
                "GaussianKernel::pdf: x is empty"
                " | kernel dim=" + std::to_string(dim_) +
                " | x dim=" + std::to_string(x.size()));

        double sq_dist = 0.0;
        for (std::size_t i = 0; i < effective_dim; ++i) {
            double d = x[i] - mu_[i];
            sq_dist += d * d;
        }
        double exponent = -0.5 * sq_dist / (sigma_ * sigma_);
        return norm_const_ * std::exp(exponent);
    }

    // Scalar overload: pdf at a single double
    double pdf(double x) const {
        return pdf(std::vector<double>{x});
    }

    // Getters (useful for debugging / pybind11 exposure)
    const std::vector<double>& mu()         const { return mu_; }
    double                     sigma()      const { return sigma_; }
    std::size_t                dim()        const { return dim_; }
    double                     norm_const() const { return norm_const_; }

private:
    std::vector<double> mu_;
    double              sigma_;
    std::size_t         dim_;
    double              norm_const_;
};