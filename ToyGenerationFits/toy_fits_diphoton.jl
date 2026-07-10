module diphoton_toy_fits

using StatsBase
using BAT
using DensityInterface
using Distributions
using Plots
using Measures
using LinearAlgebra
using Random
using LaTeXStrings

Random.seed!(42)

using ..diphoton_fit

export log_BF_results, median_log_BF_by_N_sig, all_stats_by_N_sig

# -------------------------------------------------------
# Fixed background coefficients
# -------------------------------------------------------
const p0_fixed = best.p0
const p1_fixed = best.p1
const p2_fixed = best.p2
const p3_fixed = best.p3
const p4_fixed = best.p4

# -------------------------------------------------------
# Likelihood
# -------------------------------------------------------
function log_likelihood(N_sig, toy_counts)
    x_mid  = (x_min + x_max) / 2
    x_half = (x_max - x_min) / 2
    ll = 0.0

    for i in eachindex(bin_centers)
        x = bin_centers[i]
        t = (x - x_mid) / x_half

        λ = (N_sig * s_pdf[i] +
             p0_fixed + p1_fixed*t + p2_fixed*t^2 +
             p3_fixed*t^3 + p4_fixed*t^4) * bin_width

        (!isfinite(λ) || λ <= 0) && return -Inf
        ll += logpdf(Poisson(λ), toy_counts[i])
    end

    return ll
end

function make_likelihood_for_toy(toy_counts)
    return logfuncdensity() do params
        log_likelihood(params.N_sig, toy_counts)
    end
end

# -------------------------------------------------------
# Generate toy
# -------------------------------------------------------
function generate_toy_at_N_sig(N_sig_injected)
    x_mid  = (x_min + x_max) / 2
    x_half = (x_max - x_min) / 2

    expected_counts = map(enumerate(bin_centers)) do (i, x)
        t   = (x - x_mid) / x_half
        sig = N_sig_injected * s_pdf[i]
        bkg = p0_fixed + p1_fixed*t + p2_fixed*t^2 +
              p3_fixed*t^3 + p4_fixed*t^4
        return (sig + bkg) * bin_width
    end

    return [rand(Poisson(λ)) for λ in expected_counts]
end

# -------------------------------------------------------
# Weighted polynomial fit
# -------------------------------------------------------
function weighted_polynomial_fit(x, y, yerr; degree=2)
    W = Diagonal(1.0 ./ yerr.^2)
    X = hcat([x.^i for i in 0:degree]...)
    coeffs = inv(X' * W * X) * (X' * W * y)
    return coeffs
end

# -------------------------------------------------------
# Compute log BF using parabola method
# -------------------------------------------------------
function get_log_bf(toy, N_sig_injected)

    ll = make_likelihood_for_toy(toy)

    prior = BAT.NamedTupleDist(
        N_sig = Uniform(0.0, max(1500.0, 2N_sig_injected + 1.0))
    )

    posterior = PosteriorDensity(ll, prior)

    samples = bat_sample(
        posterior,
        MCMCSampling(mcalg = MetropolisHastings(),
                     nsteps = 5*10^4,
                     nchains = 4)
    ).result

    N_vals = [s.v.N_sig for s in samples]

    h = fit(Histogram, N_vals,
            range(0.0,
                  stop=max(maximum(N_vals), 1.5*N_sig_injected + 1.0),
                  length=60))

    bw = step(h.edges[1])
    densities = h.weights ./ (length(N_vals) * bw)

    centers = (h.edges[1][1:end-1] .+ h.edges[1][2:end]) ./ 2

    nll = -log.(densities .+ eps(Float64))
    nll_err = 1.0 ./ sqrt.(max.(h.weights, 1.0))

    coeffs = weighted_polynomial_fit(centers, nll, nll_err; degree=2)
    poly(x) = coeffs[1] + coeffs[2]*x + coeffs[3]*x^2

    return max(0.0, poly(0.0) - minimum(nll))
end

# -------------------------------------------------------
# Scan
# -------------------------------------------------------
#N_sig_scan = collect(0.0:50.0:800.0)
N_sig_scan = collect(0.0:200.0:800.0)
num_toys_per_point = 15

median_log_BF_by_N_sig = Float64[]
q16_log_BF_by_N_sig = Float64[]
q84_log_BF_by_N_sig = Float64[]

all_stats_by_N_sig = Dict{Float64, Any}()

println("\n========== FINAL RESULTS ==========")

for N_inj in N_sig_scan

    results = Float64[]

    for i in 1:num_toys_per_point
        Random.seed!(i + Int(N_inj))
        toy = generate_toy_at_N_sig(N_inj)
        push!(results, get_log_bf(toy, N_inj))
    end

    clean = filter(isfinite, results)

    med = median(clean)
    q16 = quantile(clean, 0.16)
    q84 = quantile(clean, 0.84)
    avg = mean(clean)
    stdv = std(clean)

    push!(median_log_BF_by_N_sig, med)
    push!(q16_log_BF_by_N_sig, q16)
    push!(q84_log_BF_by_N_sig, q84)

    all_stats_by_N_sig[N_inj] = (
        median = med,
        mean = avg,
        std = stdv,
        q16 = q16,
        q84 = q84,
        raw = results
    )

    println("\nN_sig = $(Int(N_inj)) | median = $(round(med,digits=3))")
end

println("\n===================================\n")

# -------------------------------------------------------
# Fit (median-based)
# -------------------------------------------------------
mask = N_sig_scan .> 0

x_fit = N_sig_scan[mask]
y_fit = median_log_BF_by_N_sig[mask]

σ_fit = 0.5 .* (q84_log_BF_by_N_sig[mask] .- q16_log_BF_by_N_sig[mask])

w = 1.0 ./ (σ_fit.^2)

a = sum(w .* (x_fit.^2) .* y_fit) / sum(w .* (x_fit.^4))

fit_func(x) = a .* x.^2

println("Fitted a = ", round(a,digits=6))

T1, T2 = 4.5, 12.5

x_T1 = sqrt(T1 / a)
x_T2 = sqrt(T2 / a)

println("N_sig @4.5 = ", round(x_T1,digits=2))
println("N_sig @12.5 = ", round(x_T2,digits=2))

# -------------------------------------------------------
# Plot
# -------------------------------------------------------
lower = median_log_BF_by_N_sig .- q16_log_BF_by_N_sig
upper = q84_log_BF_by_N_sig .- median_log_BF_by_N_sig

p = scatter(
    N_sig_scan,
    median_log_BF_by_N_sig,
    yerror = (lower, upper),
    markersize = 6,
    color = :blue,
    markerstrokecolor = :black,
    xlabel = L"N_{\mathrm{sig}}",
    ylabel = L"\mathrm{median}(\log BF)",
    label = L"\mathrm{Toy\ median}"
)

xr = range(0, stop=maximum(N_sig_scan)*1.1, length=300)

plot!(
    xr, fit_func.(xr),
    lw = 3,
    color = :black,
    label = L"a N_{\mathrm{sig}}^2"
)

hline!([T1], color=:orange, linestyle=:dash, lw=2,
       label=L"\log BF=4.5")

hline!([T2], color=:green, linestyle=:dash, lw=2,
       label=L"\log BF=12.5")

display(p)

end