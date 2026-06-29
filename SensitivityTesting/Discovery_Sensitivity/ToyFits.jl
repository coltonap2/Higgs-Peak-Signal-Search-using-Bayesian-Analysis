module ToyBATFit

using CSV, DataFrames, StatsBase, BAT, Distributions, Measures, Plots, DensityInterface, LinearAlgebra
using ..HiggsAnalysis3
using ..ToyGeneration
using Random
Random.seed!(42)

export fit_toys

prior_for_toys = BAT.NamedTupleDist(
    mu = Uniform(0, 5)
)

log_likelihood_density_for_toys = logfuncdensity() do params
    log_likelihood(params.mu,
                   toys[1],
                   fit_signal,
                   fit_background)
end

# Posterior
posterior_toys = PosteriorDensity(log_likelihood_density_for_toys, prior_for_toys)

samples_toys = bat_sample(
    posterior_toys,
    MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
    # could also use 'mcalg = HamiltonianMonteCarlo'
).result

bat_report(samples_toys)
display(bat_report(samples_toys))

samples_toys_mode = mode(samples_toys)
best_fit_toys = bat_findmode(
    posterior_toys
).result

println("\nSample Mode: ", round(samples_toys_mode.mu, digits=5))
println("Best Fit Signal Strength: ", round(best_fit_toys.mu, digits=5))

mu_samples_toys = [s.mu for s in samples_toys.v]

mu_mean_toys = mean(mu_samples_toys)
mu_std_toys = std(mu_samples_toys)
mu_median_toys = median(mu_samples_toys)

mu_lower_68_t = quantile(mu_samples_toys, 0.16)
mu_upper_68_t = quantile(mu_samples_toys, 0.84)
mu_lower_95_t = quantile(mu_samples_toys, 0.025)
mu_upper_95_t = quantile(mu_samples_toys, 0.975)


# Probability Density Function
pdf_posterior_toys = histogram(mu_samples_toys,
    bins = 60,
    # normalize = :pdf,
    label = "Probability of μ given data",
    xlabel = "Signal strength (μ)",
    ylabel = "Probability Density",
    title = "Probability Density by Signal Strength (μ)",
    color = :turquoise4,
    framestyle = :box,
    size = (900, 600),
    fillalpha = 1.0,
    left_margin = 5mm,
    bottom_margin = 3mm,
    right_margin = 3mm,
    top_margin = 3mm
)

vline!([1.0], label = "SM (μ=1)", lw = 2, color = :green, linestyle = :dash)
vline!([0.0], label = "No signal (μ=0)", lw = 2, color = :red, linestyle = :dash)
vline!([mu_mean_toys], label = "Julia: μ=$(round(mu_mean_toys, digits=2))", lw = 2, linestyle = :dash, color = :yellow)

# Shade 68% credible interval
vspan!([mu_lower_68_t, mu_upper_68_t], alpha=0.2, color=:steelblue, label="68% CI")

display(pdf_posterior_toys)

prior_for_toys = BAT.NamedTupleDist(
    mu = Uniform(0, 5)
)

log_likelihood_density_for_toys = logfuncdensity() do params
    log_likelihood(params.mu,
                   toys[8],
                   fit_signal,
                   fit_background)
end

# Posterior
posterior_toys = PosteriorDensity(log_likelihood_density_for_toys, prior_for_toys)

samples_toys = bat_sample(
    posterior_toys,
    MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
    # could also use 'mcalg = HamiltonianMonteCarlo'
).result

bat_report(samples_toys)
display(bat_report(samples_toys))

samples_toys_mode = mode(samples_toys)
best_fit_toys = bat_findmode(
    posterior_toys
).result

println("\nSample Mode: ", round(samples_toys_mode.mu, digits=5))
println("Best Fit Signal Strength: ", round(best_fit_toys.mu, digits=5))

mu_samples_toys = [s.mu for s in samples_toys.v]

mu_mean_toys = mean(mu_samples_toys)
mu_std_toys = std(mu_samples_toys)
mu_median_toys = median(mu_samples_toys)

mu_lower_68_t = quantile(mu_samples_toys, 0.16)
mu_upper_68_t = quantile(mu_samples_toys, 0.84)
mu_lower_95_t = quantile(mu_samples_toys, 0.025)
mu_upper_95_t = quantile(mu_samples_toys, 0.975)


# Probability Density Function
pdf_posterior_toys = histogram(mu_samples_toys,
    bins = 60,
    # normalize = :pdf,
    label = "Probability of μ given data",
    xlabel = "Signal strength (μ)",
    ylabel = "Probability Density",
    title = "Probability Density by Signal Strength (μ)",
    color = :turquoise4,
    framestyle = :box,
    size = (900, 600),
    fillalpha = 1.0,
    left_margin = 5mm,
    bottom_margin = 3mm,
    right_margin = 3mm,
    top_margin = 3mm
)

vline!([1.0], label = "SM (μ=1)", lw = 2, color = :green, linestyle = :dash)
vline!([0.0], label = "No signal (μ=0)", lw = 2, color = :red, linestyle = :dash)
vline!([mu_mean_toys], label = "Julia: μ=$(round(mu_mean_toys, digits=2))", lw = 2, linestyle = :dash, color = :yellow)

# Shade 68% credible interval
vspan!([mu_lower_68_t, mu_upper_68_t], alpha=0.2, color=:steelblue, label="68% CI")

display(pdf_posterior_toys)


# Extract mu values from toy samples for histogram and NLL calculation
N_CE_draws = [ sample.v.mu for sample in samples_toys ] 
# Building a histogram of the mu values to estimate the pdf and then compute the NLL
lowerb = minimum(N_CE_draws)
upperb = maximum(N_CE_draws)
edges = range(lowerb, stop=upperb, length=50)
hist_obj   = fit(Histogram, N_CE_draws, edges)
counts     = hist_obj.weights # how many toy experiments fall into each bin
# converts histogram counts to pdf
bin_centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
total      = length(N_CE_draws)
bw         = step(edges)
densities  = counts ./ (total * bw)
# Calculate Negative Log-Likelihood (NLL) from the densities
nll        = -log.(densities .+ eps(Float64))
# saftey to avoid log(0)
eps_val = eps(Float64)
# relative error bars for NLL based on Poisson statistics of the counts in each bin
nll_errors = 1.0 ./ sqrt.(max.(counts, eps_val))
nll_std = std(nll)

inverted_plot = scatter(
    bin_centers,
    nll,
    ylim = (0, 15),
    yerror = nll_errors,  # Error bars from Poisson statistics
    xlabel = "mu (Signal Strength)",
    ylabel = "Negative Log-Likelihood (NLL)",
    title = "NLL vs mu with Toy Fits",
    framestyle = :box,
    size = (900, 600),
    left_margin = 5mm,
    bottom_margin = 3mm,
    right_margin = 3mm,
    top_margin = 3mm
)

display(inverted_plot)

# Function tries polynomial degrees 4, 6, 8
# Fits using weighted least squares with weights from yerr
# Rejects bad fits based on local maxima, bad chi squared/p-value, or too many parameters hitting upper bound

function weighted_polynomial_fit(x, y, yerr; degrees=[2,4,6], x0=0.0, upper_bound=samples_toys_mode.mu + 1 * nll_std)
    best_p = -1.0 # Initial placeholder p value
    best_result = nothing
    W = Diagonal(1.0 ./ yerr.^2) # Weights for the fit
    for deg in degrees
        xi = x .- x0 # Shift x values to improve numerical stability
        X = hcat([xi.^i for i in 0:deg]...) # Design matrix for polynomial terms
        cov = inv(X' * W * X) # Covariance matrix of the fit parameters
        coeffs = cov * (X' * W * y) # Fit coefficients
        y_pred = X * coeffs # Fitted values
        chi2 = sum(((y - y_pred) ./ yerr).^2) # Chi-squared statistic
        dof = length(y) - (deg + 1) # Degrees of freedom
        p_value = 1 - cdf(Chisq(dof), chi2) # p-value for goodness of fit
        # Reject bad statistical fits
        if p_value < 0.001
            continue
        end
        # Reject fits with local maxima
        deriv_coeffs = [i * coeffs[i+1] for i in 1:deg] # Coefficients of the derivative
        has_max = false # placeholder for local max check
        xs = range(0, stop=upper_bound, length=100)
            for x_test in xs
                xt = x_test - x0
            # derivative at x
            d = sum(deriv_coeffs[i] * xt^(i-1) for i in 1:length(deriv_coeffs))

            # check sign change
            xt_before = (x_test - 0.01) - x0
            xt_after  = (x_test + 0.01) - x0

            d_before = sum(deriv_coeffs[i] * xt_before^(i-1) for i in 1:length(deriv_coeffs))
            d_after  = sum(deriv_coeffs[i] * xt_after^(i-1) for i in 1:length(deriv_coeffs))

            if d_before > 0 && d_after < 0
                has_max = true
                break
            end
        end
        if has_max
            continue
        end
        if p_value > best_p
            best_p = p_value
            best_result = (coeffs=coeffs, chi2=chi2, dof=dof, p_value=p_value, degree=deg)
        end
    end

    return best_result
end

weighted_polynomial_fit_result = weighted_polynomial_fit(bin_centers, nll, nll_errors)
if weighted_polynomial_fit_result !== nothing
    println("Best Fit Degree: ", weighted_polynomial_fit_result.degree)
    println("Chi-squared: ", weighted_polynomial_fit_result.chi2)
    println("Degrees of Freedom: ", weighted_polynomial_fit_result.dof)
    println("Reduced Chi-squared: ", weighted_polynomial_fit_result.chi2 / weighted_polynomial_fit_result.dof)
    println("p-value: ", weighted_polynomial_fit_result.p_value)
    println("Coefficients: ", weighted_polynomial_fit_result.coeffs)
else
    println("No acceptable fit found.")
end

polynomial_fit = x -> sum(weighted_polynomial_fit_result.coeffs[i] * (x)^(i-1) for i in 1:length(weighted_polynomial_fit_result.coeffs))

# Explicitly display the fit function
println("\n\n=== POLYNOMIAL FIT FUNCTION ===")
println("Degree: ", weighted_polynomial_fit_result.degree)
print("Function: f(x) = ")
for i in 1:length(weighted_polynomial_fit_result.coeffs)
    coeff = weighted_polynomial_fit_result.coeffs[i]
    power = i - 1
    if power == 0
        print("$(round(coeff, digits=6))")
    elseif power == 1
        print("+ $(round(coeff, digits=6))*(x - $(round(samples_toys_mode.mu, digits=4)))")
    else
        print("+ $(round(coeff, digits=6))*(x - $(round(samples_toys_mode.mu, digits=4)))^$power")
    end
    if i < length(weighted_polynomial_fit_result.coeffs)
        print(" ")
    end
end
println("\n\nCoefficients: ")
for i in 1:length(weighted_polynomial_fit_result.coeffs)
    println("  c[$i] = $(weighted_polynomial_fit_result.coeffs[i])")
end
println("================================\n")

# Plot the fit function alongside the data
x_range = range(0, stop=samples_toys_mode.mu, length=200)
y_fit = polynomial_fit.(x_range)

# Reduce range to focus on region around mode
filter_indices = (bin_centers .>= 0) .& (bin_centers .<= samples_toys_mode.mu)
bin_centers_filtered = bin_centers[filter_indices]
nll_filtered = nll[filter_indices]
nll_errors_filtered = nll_errors[filter_indices]

plot_with_fit = scatter(
    bin_centers_filtered,
    nll_filtered,
    yerror = nll_errors_filtered,
    label = "Data",
    xlabel = "mu (Signal Strength)",
    ylabel = "Negative Log-Likelihood (NLL)",
    title = "NLL vs mu with Polynomial Fit",
    framestyle = :box,
    size = (900, 600),
    left_margin = 5mm,
    bottom_margin = 3mm,
    right_margin = 3mm,
    top_margin = 3mm,
    legend = :topright
)

plot!(plot_with_fit, x_range, y_fit,
    label = "Polynomial Fit (degree $(weighted_polynomial_fit_result.degree))",
    ylim = (0, 15),
    lw = 3,
    color = :red,
    linestyle = :solid
)

display(plot_with_fit)

H_0 = -polynomial_fit(0.0)
H_1 = -minimum(nll)

print("Log BF = H_1 - H_0 = ", round(H_1 - H_0, digits=4), "\n")

function get_log_bf(current_toy)
    # 1. Define Likelihood for THIS toy
    log_likelihood_density = logfuncdensity() do params
        log_likelihood(params.mu, current_toy, fit_signal, fit_background)
    end
    posterior = PosteriorDensity(log_likelihood_density, prior_for_toys)

    # 2. Sample
    samples = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^4, nchains = 2)).result
    
    # 3. Process samples to NLL
    mu_vals = [s.v.mu for s in samples]
    h = fit(Histogram, mu_vals, range(minimum(mu_vals), stop=maximum(mu_vals), length=40))
    
    # PDF to NLL conversion
    bw = step(h.edges[1])
    densities = h.weights ./ (length(mu_vals) * bw)
    bin_centers = (h.edges[1][1:end-1] .+ h.edges[1][2:end]) ./ 2
    nll = -log.(densities .+ eps(Float64))
    nll_err = 1.0 ./ sqrt.(max.(h.weights, eps(Float64)))

    # 4. Fit Polynomial
    fit_res = weighted_polynomial_fit(bin_centers, nll, nll_err; degrees=[2,4,6])
    
    if fit_res !== nothing
        poly = x -> sum(fit_res.coeffs[i] * x^(i-1) for i in 1:length(fit_res.coeffs))
        
        # Calculate Log Bayes Factor (H1 - H0)
        # Note: In NLL terms, H = -NLL
        h0 = -poly(0.0)
        h1 = -minimum(nll)
        return max(0, h1 - h0) # Ensure non-negative Log BF
    else
        return NaN # Return Not-a-Number if fit fails
    end
end

# --- RUN THE LOOP ---
num_toys = 100
log_BF_results = []

println("Starting fits for $num_toys toys...")
for i in 1:num_toys
    Random.seed!(i)
    lbf = get_log_bf(toys[i])
    push!(log_BF_results, lbf)
    println("Toy $i: Log BF = ", round(lbf, digits=4))
end

println("\nFinal Results: ", log_BF_results)

# 1. Define positions where you want the ticks to appear
tick_positions = [collect(0:2:18)..., 20.5]

# 2. Define the corresponding labels as strings
tick_labels = [string.(0:2:18)..., "20+"]

# After running your loop for 100 toys
final_graph = histogram(log_BF_results, 
    bins = 0:1:21, 
    xticks = (tick_positions, tick_labels),
    xlabel = "Log Bayes Factor",
    title = "Evidence Distribution with Overflow Bin",
    label = "Toy Fits",
    color = :steelblue,
    fillalpha = 0.8)

display(final_graph)

# Filter out NaNs and then take the mean
clean_results = filter(x -> !isnan(x) && x != 0.0, log_BF_results)
avg_log_bf = mean(clean_results)

println("Mean Log BF (excluding NaNs): ", avg_log_bf)
end