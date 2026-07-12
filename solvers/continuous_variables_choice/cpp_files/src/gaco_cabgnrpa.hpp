#pragma once

#include <vector>
#include <functional>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <memory>

// pagmo headers
#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

#include "cnrpa.hpp"          // Policy, RNG, adapt_policy, helpers
#include "gaussian_kernel.hpp"
// pagmo_udp.hpp included via bindings.cpp — not needed here

// ---------------------------------------------------------------------------
// sample_mixture_1d
// Draws one sample from a 2-component 1D Gaussian mixture.
// weight1 = probability of drawing from component 1.
// ---------------------------------------------------------------------------
inline double sample_mixture_1d(
    double mu1, double sigma1,
    double mu2, double sigma2,
    double weight1,
    RNG&   rng)
{
    bool use1 = (rng.uniform(0.0, 1.0) < weight1);
    return use1 ? rng.normal(mu1, sigma1) : rng.normal(mu2, sigma2);
}

// ---------------------------------------------------------------------------
// fit_gaussian_from_density
// Weighted mean + std from (x, density) pairs.
// std_factor multiplies the computed sigma (mirrors the added Python argument).
// Returns (mu, sigma).
// ---------------------------------------------------------------------------
inline std::pair<double, double> fit_gaussian_from_density(
    const std::vector<double>& x,
    const std::vector<double>& density,
    double                     std_factor)
{
    if (x.size() != density.size() || x.empty())
        throw std::invalid_argument("fit_gaussian_from_density: size mismatch");

    double sum_d = 0.0;
    for (double d : density) sum_d += d;

    std::vector<double> p(density.size());
    for (std::size_t i = 0; i < density.size(); ++i)
        p[i] = density[i] / sum_d;

    double mu = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
        mu += p[i] * x[i];

    double var = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        double d = x[i] - mu;
        var += p[i] * d * d;
    }

    return {mu, std_factor * std::sqrt(var)};
}

// ---------------------------------------------------------------------------
// BiasState: per-dimension bias mean and std derived from GACO population
// ---------------------------------------------------------------------------
struct BiasState {
    std::vector<double> bias;
    std::vector<double> bias_std;
};

// ---------------------------------------------------------------------------
// GACOState
// Owns the pagmo algorithm and population directly in C++.
// No Python callbacks needed — GACO runs purely in C++.
// The population is copyable (pagmo::population has value semantics),
// so save/restore checkpoints are just copies.
// ---------------------------------------------------------------------------
struct GACOState {
    pagmo::algorithm  algo;
    pagmo::population pop;
    pagmo::population checkpoint;   // saved snapshot for restore

    GACOState(pagmo::algorithm a, pagmo::population p)
        : algo(std::move(a))
        , pop(std::move(p))
        , checkpoint(pop)   // initialize checkpoint to initial pop
    {}

    void evolve() {
        pop = algo.evolve(pop);
    }

    void save() {
        checkpoint = pop;   // pagmo::population is copyable
    }

    void restore() {
        pop = checkpoint;
    }

    const std::vector<pagmo::vector_double>& get_x() const {
        return pop.get_x();
    }

    std::vector<double> get_f_flat() const {
        const auto& f = pop.get_f();
        std::vector<double> flat;
        flat.reserve(f.size());
        for (const auto& row : f)
            flat.push_back(row.empty() ? 0.0 : row[0]);
        return flat;
    }
};

// ---------------------------------------------------------------------------
// gaco_cabgnrpa_playout
// Mirrors Python's gaco_cabgnrpa_playout().
// ---------------------------------------------------------------------------
inline std::pair<std::vector<double>, std::vector<double>>
gaco_cabgnrpa_playout(
    const Policy&    policy,
    const BiasState& bias_state,
    const Bounds&    bounds,
    double           std_factor,
    double           tau,
    double           gaussian_kernel_threshold,
    RNG&             rng)
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
            if (has_policy && policy[0].is_scalar) {
                double policy_sample = denormalize(
                    rng.normal(policy[0].scalar_value, std_factor), lo, hi);

                chosen_value = truncate(
                    sample_mixture_1d(
                        policy_sample,
                        std_factor * (hi - lo),
                        bias_state.bias[advancement],
                        bias_state.bias_std[advancement],
                        1.0 / tau,
                        rng),
                    lo, hi);
            } else {
                chosen_value = rng.uniform(lo, hi);
            }
            states_sequence.push_back(normalize(chosen_value, lo, hi));

        } else {
            double policy_candidate = 0.0;
            bool   has_candidate    = false;

            if (has_policy && advancement < policy.size()
                && !policy[advancement].table.empty())
            {
                // 1D kernel centered at first element (matches numpy truncation)
                GaussianKernel kernel(
                    std::vector<double>{states_sequence[0]}, std_factor);

                std::vector<double> weights, norm_values;
                for (const auto& [key, entry] : policy[advancement].table) {
                    double w = kernel.pdf(std::vector<double>{entry.norm_value});
                    if (w >= gaussian_kernel_threshold) {
                        weights.push_back(w);
                        norm_values.push_back(entry.norm_value);
                    }
                }

                if (!norm_values.empty()) {
                    double sum_w = 0.0;
                    for (double w : weights) sum_w += w;
                    for (double& w : weights) w /= sum_w;

                    double wm = 0.0;
                    for (std::size_t i = 0; i < norm_values.size(); ++i)
                        wm += weights[i] * norm_values[i];

                    policy_candidate = denormalize(
                        rng.normal(wm, std_factor), lo, hi);
                    has_candidate = true;
                }
            }

            if (!has_candidate)
                policy_candidate = rng.uniform(lo, hi);

            chosen_value = truncate(
                sample_mixture_1d(
                    policy_candidate,
                    std_factor * (hi - lo) / 2.0,
                    bias_state.bias[advancement],
                    bias_state.bias_std[advancement],
                    1.0 / tau,
                    rng),
                lo, hi);

            states_sequence.push_back(normalize(chosen_value, lo, hi));
        }

        values_sequence.push_back(chosen_value);
    }

    return {values_sequence, states_sequence};
}

// ---------------------------------------------------------------------------
// GACOCNRPAResult
// ---------------------------------------------------------------------------
struct GACOCNRPAResult {
    std::vector<double> best_values_sequence;
    std::vector<double> best_states_sequence;
    double              best_value;
    std::vector<double> best_values_list;
    std::vector<double> time_list;
};

// ---------------------------------------------------------------------------
// run_gaco_cabgnrpa  (recursive)
// ---------------------------------------------------------------------------
inline GACOCNRPAResult run_gaco_cabgnrpa(
    const Evaluator&                       evaluator,
    Policy                                 policy,
    GACOState&                             gaco,
    int                                    level,
    int                                    n_policies,
    double                                 zeta,
    const Bounds&                          bounds,
    int                                    current_iteration,
    double                                 learning_rate,
    double                                 timeout,
    double                                 unfeasibility_value,
    double                                 gaussian_kernel_threshold,
    double                                 tau,
    RNG&                                   rng,
    std::chrono::steady_clock::time_point  start_time,
    std::vector<double>                    best_values_sequence = {},
    std::vector<double>                    best_states_sequence = {},
    double                                 best_value =
                                               std::numeric_limits<double>::max(),
    std::vector<double>                    best_values_list = {},
    std::vector<double>                    time_list        = {})
{
    auto elapsed = [&]() -> double {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    };

    if (level == 0) {
        // --- Level 0: evolve GACO, compute bias, playout ---

        // GACO evolve: pure C++ pagmo, no GIL needed
        gaco.evolve();

        // Extract population
        const auto& bias_x  = gaco.get_x();
        auto        bias_f  = gaco.get_f_flat();

        // Density weights = 1 / fitness, normalized
        std::vector<double> density(bias_f.size());
        for (std::size_t i = 0; i < bias_f.size(); ++i)
            density[i] = (bias_f[i] > 0.0) ? 1.0 / bias_f[i] : 0.0;

        double sum_d = 0.0;
        for (double d : density) sum_d += d;
        if (sum_d > 0.0)
            for (double& d : density) d /= sum_d;

        // Fit per-dimension Gaussian
        BiasState bias_state;
        bias_state.bias.resize(bounds.size());
        bias_state.bias_std.resize(bounds.size());

        for (std::size_t dim = 0; dim < bounds.size(); ++dim) {
            std::vector<double> x_dim(bias_x.size());
            for (std::size_t i = 0; i < bias_x.size(); ++i)
                x_dim[i] = bias_x[i][dim];

            auto [mu, sigma] = fit_gaussian_from_density(x_dim, density, zeta);
            bias_state.bias[dim]     = mu;
            bias_state.bias_std[dim] = sigma;
        }

        // Playout
        double std_factor = 0.01 + 1.0 /
            std::sqrt(static_cast<double>(current_iteration) + 1.0);

        auto [values_sequence, states_sequence] = gaco_cabgnrpa_playout(
            policy, bias_state, bounds, std_factor, tau,
            gaussian_kernel_threshold, rng);

        double fitness = evaluator(values_sequence);

        return GACOCNRPAResult{
            values_sequence, states_sequence, fitness,
            best_values_list, time_list};

    } else {
        // --- Recursive case ---
        // GACO population is NOT reset between iterations — it accumulates
        // across all n_policies iterations, matching Python's behavior where
        // deepcopy(biases_values) copies once before the loop and is never
        // reset during it.
        Policy current_policy = policy;

        for (int iter = 0; iter < n_policies; ++iter) {
            GACOCNRPAResult sub = run_gaco_cabgnrpa(
                evaluator, current_policy, gaco,
                level - 1, n_policies, zeta, bounds,
                iter, learning_rate, timeout,
                unfeasibility_value, gaussian_kernel_threshold, tau,
                rng, start_time,
                best_values_sequence, best_states_sequence,
                best_value, best_values_list, time_list);

            if (sub.best_value < best_value) {
                best_value           = sub.best_value;
                best_values_sequence = sub.best_values_sequence;
                best_states_sequence = sub.best_states_sequence;
            }

            best_values_list = std::move(sub.best_values_list);
            time_list        = std::move(sub.time_list);

            double t = elapsed();

            if (best_value < unfeasibility_value) {
                adapt_policy(
                    best_values_sequence, best_states_sequence,
                    current_policy, learning_rate, bounds,
                    gaussian_kernel_threshold);

                best_values_list.push_back(best_value);
                time_list.push_back(t);
            }

            if (t > timeout) break;

            // NOTE: do NOT restore GACO population here.
            // Python's deepcopy(biases_values) copies the population ONCE
            // before the loop but never resets it during the loop —
            // GACO accumulates improvements across all n_policies iterations
            // at this level. We match that by never calling gaco.restore()
            // inside the loop.
        }

        // Use best_value directly — already the best fitness found,
        // no need to re-evaluate (mirrors corrected Python behavior).
        return GACOCNRPAResult{
            best_values_sequence, best_states_sequence, best_value,
            best_values_list, time_list};
    }
}

// ---------------------------------------------------------------------------
// gaco_cabgnrpa  — top-level entry point
// Constructs pagmo::algorithm and pagmo::population entirely in C++.
// Only the fitness() callback still touches Python (unavoidable since
// the trajectory evaluator lives in pykep Python).
// ---------------------------------------------------------------------------
inline GACOCNRPAResult gaco_cabgnrpa(
    const Evaluator&            evaluator,
    pagmo::problem&             prob,           // pre-built with GIL held in bindings.cpp
    const pagmo::vector_double& lb,
    const pagmo::vector_double& ub,
    int                         level,
    int                         n_policies,
    double                      zeta,
    const Bounds&               bounds,
    double                      learning_rate,
    double                      timeout,
    double                      tau,
    double                      unfeasibility_value,
    double                      gaussian_kernel_threshold,
    unsigned                    kernel_size,
    unsigned                    n_generations,
    double                      elitism_factor,
    unsigned                    random_seed)
{
    // prob is pre-built in bindings.cpp with GIL held — do not reconstruct here.

    // Build pagmo GACO algorithm
    pagmo::algorithm algo{pagmo::gaco(
        n_generations,
        kernel_size,
        elitism_factor,   // q
        1e9,              // oracle (large = unknown best)
        0.01,             // acc
        1u,               // threshold
        7u,               // n_gen_mark
        0.0,              // focus
        true,             // memory = true
        random_seed)};

    // Build initial population
    pagmo::population pop{prob, kernel_size, random_seed};

    // Construct GACOState
    GACOState gaco_state(std::move(algo), std::move(pop));

    RNG  rng(random_seed);
    auto start = std::chrono::steady_clock::now();

    return run_gaco_cabgnrpa(
        evaluator, Policy{}, gaco_state,
        level, n_policies, zeta, bounds,
        0, learning_rate, timeout,
        unfeasibility_value, gaussian_kernel_threshold, tau,
        rng, start);
}