#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <functional>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <limits>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "gaussian_kernel.hpp"

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------
using Bound    = std::pair<double, double>;
using Bounds   = std::vector<Bound>;
using Evaluator = std::function<double(const std::vector<double>&)>;

// ---------------------------------------------------------------------------
// PolicyEntry
//
// At advancement == 0: scalar mean for departure epoch (is_scalar = true)
// At advancement >= 1: table mapping state-key -> (normalized_value, state_vector)
//
// We store the actual state vector alongside each entry so we can correctly
// evaluate GaussianKernel.pdf(key) as in the Python version, where keys ARE
// the state tuples.
// ---------------------------------------------------------------------------
struct PolicyTableEntry {
    double              norm_value;    // normalized decision value
    std::vector<double> state_vec;     // the state prefix that produced this key
};

struct PolicyEntry {
    bool   is_scalar    = false;
    double scalar_value = 0.5;

    // key -> {norm_value, state_vec}
    // key = code(state_vec), used for fast lookup
    std::unordered_map<std::string, PolicyTableEntry> table;
};

using Policy = std::vector<PolicyEntry>;

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------
inline double normalize(double value, double lo, double hi) {
    return (value - lo) / (hi - lo);
}

inline double denormalize(double norm, double lo, double hi) {
    return norm * (hi - lo) + lo;
}

inline double truncate(double value, double lo, double hi) {
    return std::max(lo, std::min(hi, value));
}

// Mirrors Python's code():
//   rounds each value and returns a tuple-string like "(0.25, 0.73)"
// Used both as dict key and as the state representation passed to pdf().
inline std::string code(const std::vector<double>& states) {
    std::ostringstream oss;
    oss << "(";
    for (std::size_t i = 0; i < states.size(); ++i) {
        // Round to 2 decimal places (matching Python's round() on normalized vals)
        double rounded = std::round(states[i] * 100.0) / 100.0;
        oss << std::fixed << std::setprecision(2) << rounded;
        if (i + 1 < states.size()) oss << ", ";
    }
    oss << ")";
    return oss.str();
}

// Parse a code() string back to a vector<double>
// "(0.25, 0.73)" -> {0.25, 0.73}
inline std::vector<double> decode(const std::string& key) {
    std::vector<double> result;
    std::string s = key;
    // Remove parentheses
    s.erase(std::remove(s.begin(), s.end(), '('), s.end());
    s.erase(std::remove(s.begin(), s.end(), ')'), s.end());
    std::istringstream iss(s);
    std::string token;
    while (std::getline(iss, token, ',')) {
        // Trim whitespace
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);
        if (!token.empty())
            result.push_back(std::stod(token));
    }
    return result;
}

// ---------------------------------------------------------------------------
// RNG: mirrors numpy.random.default_rng
// ---------------------------------------------------------------------------
struct RNG {
    std::mt19937_64                        engine;
    std::uniform_real_distribution<double> uniform_dist{0.0, 1.0};
    std::normal_distribution<double>       normal_dist{0.0, 1.0};

    explicit RNG(uint64_t seed = 0)
        : engine(seed == 0 ? std::random_device{}() : seed) {}

    double uniform(double lo, double hi) {
        return lo + (hi - lo) * uniform_dist(engine);
    }

    // mirrors: np.random.normal(mean, std)
    double normal(double mean, double std) {
        return mean + std * normal_dist(engine);
    }
};

// ---------------------------------------------------------------------------
// policy_playout
// Mirrors Python's policy_playout() exactly.
// ---------------------------------------------------------------------------
inline std::pair<std::vector<double>, std::vector<double>>
policy_playout(
    const Policy& policy,
    const Bounds& bounds,
    double        std_factor,
    double        gaussian_kernel_threshold,
    RNG&          rng)
{
    const bool has_policy = !policy.empty();

    std::vector<double> values_sequence;
    std::vector<double> states_sequence;
    values_sequence.reserve(bounds.size());
    states_sequence.reserve(bounds.size());

    for (std::size_t advancement = 0; advancement < bounds.size(); ++advancement) {
        const double lo = bounds[advancement].first;
        const double hi = bounds[advancement].second;
        double chosen_value = 0.0;

        if (advancement == 0) {
            // --- Departure epoch ---
            if (has_policy && policy[0].is_scalar) {
                double norm_sample = rng.normal(policy[0].scalar_value, std_factor);
                // Python: round(denormalize(...))
                chosen_value = std::round(denormalize(norm_sample, lo, hi));
            } else {
                chosen_value = rng.uniform(lo, hi);
            }
            states_sequence.push_back(normalize(chosen_value, lo, hi));

        } else {
            if (has_policy && advancement < policy.size()
                && !policy[advancement].table.empty())
            {
                const auto& table = policy[advancement].table;

                // In Python:
                //   gaussian_kernel = GaussianKernel(states_sequence, sigma=std_factor)
                //   value = current_policy[key]   <- a scalar (float/int)
                //   weight = gaussian_kernel.pdf(value)
                //     -> pdf wraps scalar as [value], giving a 1-element vector
                //     -> but kernel center is states_sequence (multi-dim)
                //     -> numpy silently computes dot product on min(dim) elements
                //        i.e. only the FIRST element of states_sequence matters
                //
                // To replicate this faithfully and avoid the dimension mismatch,
                // we instantiate the kernel as 1D, centered at the FIRST element
                // of states_sequence (the only dimension numpy actually uses).
                // This matches what Python computes without relying on numpy's
                // silent broadcasting behavior.
                GaussianKernel kernel(
                    std::vector<double>{states_sequence[0]}, std_factor);

                std::vector<double> weights;
                std::vector<double> norm_values;

                for (const auto& [key, entry] : table) {
                    // pdf called with the stored scalar norm_value,
                    // wrapped as a 1-element vector — matches Python's
                    // pdf([scalar]) behavior after the isinstance() wrapping.
                    double w = kernel.pdf(std::vector<double>{entry.norm_value});
                    if (w >= gaussian_kernel_threshold) {
                        weights.push_back(w);
                        norm_values.push_back(entry.norm_value);
                    }
                }

                if (!norm_values.empty()) {
                    // Normalize weights
                    double sum_w = 0.0;
                    for (double w : weights) sum_w += w;
                    for (double& w : weights) w /= sum_w;

                    // Weighted mean
                    double weighted_mean = 0.0;
                    for (std::size_t i = 0; i < norm_values.size(); ++i)
                        weighted_mean += weights[i] * norm_values[i];

                    // Sample and truncate
                    double norm_sample = rng.normal(weighted_mean, std_factor);
                    chosen_value = truncate(
                        denormalize(norm_sample, lo, hi),
                        lo, hi);
                } else {
                    chosen_value = rng.uniform(lo, hi);
                }
            } else {
                chosen_value = rng.uniform(lo, hi);
            }

            states_sequence.push_back(normalize(chosen_value, lo, hi));
        }

        values_sequence.push_back(chosen_value);
    }

    return {values_sequence, states_sequence};
}

// ---------------------------------------------------------------------------
// adapt_policy
// Mirrors Python's adapt_policy() exactly, including the GaussianKernel
// centered at current_key evaluated at each other key in the table.
// ---------------------------------------------------------------------------
inline void adapt_policy(
    const std::vector<double>& best_values_sequence,
    const std::vector<double>& best_states_sequence,
    Policy&                    policy,
    double                     learning_rate,
    const Bounds&              bounds,
    double                     gaussian_kernel_threshold)
{
    const bool has_policy = !policy.empty();

    for (std::size_t advancement = 0;
         advancement < best_values_sequence.size(); ++advancement)
    {
        const double lo      = bounds[advancement].first;
        const double hi      = bounds[advancement].second;
        const double element = best_values_sequence[advancement];
        const double norm_el = normalize(element, lo, hi);

        if (!has_policy) {
            // --- Cold init ---
            if (advancement >= policy.size())
                policy.resize(advancement + 1);

            if (advancement == 0) {
                policy[0].is_scalar    = true;
                policy[0].scalar_value = norm_el;
            } else {
                std::vector<double> prefix(
                    best_states_sequence.begin(),
                    best_states_sequence.begin() + advancement);
                std::string key = code(prefix);
                policy[advancement].is_scalar = false;
                policy[advancement].table[key] = {norm_el, prefix};
            }

        } else {
            if (advancement >= policy.size())
                policy.resize(advancement + 1);

            if (advancement == 0) {
                // Gradient step on scalar mean
                policy[0].scalar_value +=
                    learning_rate * (norm_el - policy[0].scalar_value);
            } else {
                std::vector<double> prefix(
                    best_states_sequence.begin(),
                    best_states_sequence.begin() + advancement);
                std::string current_key = code(prefix);
                auto& table = policy[advancement].table;

                // Update or insert current key
                if (table.count(current_key)) {
                    double prev = table[current_key].norm_value;
                    table[current_key].norm_value =
                        prev + learning_rate * (norm_el - prev);
                } else {
                    table[current_key] = {norm_el, prefix};
                }

                // Propagate to nearby keys via Gaussian kernel.
                // Mirrors Python:
                //   gaussian_kernel = GaussianKernel(current_key, 0.2)
                //   for key in policy[advancement].keys():
                //       weight = gaussian_kernel.pdf(key)
                //
                // Here current_key encodes `prefix` (a state vector).
                // The kernel is centered at that state vector.
                // pdf(key) evaluates at the KEY's state vector.
                GaussianKernel kernel(prefix, 0.2);

                for (auto& [key, entry] : table) {
                    if (key == current_key) continue;

                    // Guard: only call pdf if dimensions match the kernel center.
                    // Entries inserted at a different recursion depth or from a
                    // different branch may have a state_vec of a different size.
                    if (entry.state_vec.size() != prefix.size()) continue;

                    double w = kernel.pdf(entry.state_vec);
                    if (w >= gaussian_kernel_threshold) {
                        entry.norm_value +=
                            learning_rate * w * (norm_el - entry.norm_value);
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// CNRPAResult
// ---------------------------------------------------------------------------
struct CNRPAResult {
    std::vector<double> best_values_sequence;
    std::vector<double> best_states_sequence;
    double              best_value;
    std::vector<double> best_values_list;
    std::vector<double> time_list;
};

// ---------------------------------------------------------------------------
// run_cnrpa  (recursive)
// Mirrors Python's run_cnrpa() exactly.
// ---------------------------------------------------------------------------
inline CNRPAResult run_cnrpa(
    const Evaluator&                       evaluator,
    Policy                                 policy,
    int                                    level,
    int                                    n_policies,
    const Bounds&                          bounds,
    int                                    current_iteration,
    double                                 learning_rate,
    double                                 timeout,
    double                                 unfeasibility_value,
    double                                 gaussian_kernel_threshold,
    RNG&                                   rng,
    std::chrono::steady_clock::time_point  start_time,
    std::vector<double>                    best_values_sequence = {},
    std::vector<double>                    best_states_sequence = {},
    double                                 best_value           =
                                               std::numeric_limits<double>::max(),
    std::vector<double>                    best_values_list     = {},
    std::vector<double>                    time_list            = {})
{
    auto elapsed = [&]() -> double {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    };

    if (level == 0) {
        // --- Base case: single playout ---
        // std_factor mirrors: 0.01 + 1/sqrt(current_iteration + 2)
        double std_factor = 0.01 + 1.0 /
            std::sqrt(static_cast<double>(current_iteration) + 2.0);

        auto [values_sequence, states_sequence] = policy_playout(
            policy, bounds, std_factor, gaussian_kernel_threshold, rng);

        double fitness = evaluator(values_sequence);

        return CNRPAResult{
            values_sequence,
            states_sequence,
            fitness,
            best_values_list,
            time_list
        };

    } else {
        // --- Recursive case ---
        Policy current_policy = policy;   // C++ value copy — no deepcopy needed

        for (int iter = 0; iter < n_policies; ++iter) {
            CNRPAResult sub = run_cnrpa(
                evaluator,
                current_policy,
                level - 1,
                n_policies,
                bounds,
                iter,                     // current_iteration at sub-level
                learning_rate,
                timeout,
                unfeasibility_value,
                gaussian_kernel_threshold,
                rng,
                start_time,
                best_values_sequence,
                best_states_sequence,
                best_value,
                best_values_list,
                time_list);

            if (sub.best_value < best_value) {
                best_value           = sub.best_value;
                best_values_sequence = sub.best_values_sequence;
                best_states_sequence = sub.best_states_sequence;
            }

            // Carry accumulated lists forward
            best_values_list = std::move(sub.best_values_list);
            time_list        = std::move(sub.time_list);

            double t = elapsed();

            if (best_value < unfeasibility_value) {
                adapt_policy(
                    best_values_sequence,
                    best_states_sequence,
                    current_policy,
                    learning_rate,
                    bounds,
                    gaussian_kernel_threshold);

                best_values_list.push_back(best_value);
                time_list.push_back(t);
            }

            if (t > timeout) break;
        }

        return CNRPAResult{
            best_values_sequence,
            best_states_sequence,
            best_value,
            best_values_list,
            time_list
        };
    }
}

// ---------------------------------------------------------------------------
// cnrpa  — top-level entry point
// ---------------------------------------------------------------------------
inline CNRPAResult cnrpa(
    const Evaluator& evaluator,
    int              level,
    int              n_policies,
    const Bounds&    bounds,
    double           learning_rate,
    double           timeout,
    double           unfeasibility_value,
    double           gaussian_kernel_threshold,
    uint64_t         seed = 0)
{
    RNG  rng(seed);
    auto start = std::chrono::steady_clock::now();

    return run_cnrpa(
        evaluator,
        Policy{},
        level,
        n_policies,
        bounds,
        0,
        learning_rate,
        timeout,
        unfeasibility_value,
        gaussian_kernel_threshold,
        rng,
        start);
}