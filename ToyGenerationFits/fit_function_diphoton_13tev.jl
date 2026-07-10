module diphoton_fit

using CSV
using DataFrames
using StatsBase
using SpecialFunctions
using BAT
using DensityInterface
using Distributions
using IntervalSets
using Plots
using Optim
using DelimitedFiles

export fit_function, 
       best, 
       s_pdf, 
       bin_centers, 
       bin_width, 
       x_min, 
       x_max, 
       hist_data, 
       samples,
       posterior
       

data_myy = vec(readdlm("data_myy.dat", comments=true))
mc = readdlm("mc_signal_myy_weights.dat", comments=true)
mc_myy = mc[:, 1]
mc_w   = mc[:, 2]

x_min, x_max, bin_width = 100.0, 160.0, 2.0
edges = x_min:bin_width:x_max

hist_data = append!(Histogram(edges), data_myy)
n = hist_data.weights
bin_centers = (hist_data.edges[1][1:end-1] .+ hist_data.edges[1][2:end]) ./ 2

s = zeros(length(bin_centers))
for k in eachindex(mc_myy)
    idx = searchsortedlast(collect(edges), mc_myy[k])
    if 1 <= idx <= length(s)
        s[idx] += mc_w[k]
    end
end

s_pdf = s ./ (sum(s) * bin_width)

x_mid = (x_min + x_max) / 2
x_half = (x_max - x_min) / 2

function fit_function(par::NamedTuple{(:N_sig, :p0, :p1, :p2, :p3, :p4)}, x::Real)
    x_min <= x <= x_max || return 0.0
    t = (x - x_mid) / x_half
    #idx = clamp(Int(floor((x - x_min) / bin_width)) + 1, 1, length(s_pdf))
    #sig = par.N_sig * s_pdf[idx]
    sig = par.N_sig * pdf(Normal(124.96, 4.56), x)
    bkg = par.p0 + par.p1*t + par.p2*t^2 + par.p3*t^3 + par.p4*t^4
    return (sig + bkg) * bin_width
end

# function fit_function(par, x::Real)
#     x_min <= x <= x_max || return 0.0
#     t = (x - x_mid) / x_half
#     idx = clamp(Int(floor((x - x_min) / bin_width)) + 1, 1, length(s_pdf))
#     sig = par.N_sig * s_pdf[idx]
#     bkg = par.p0 + par.p1*t + par.p2*t^2 + par.p3*t^3 + par.p4*t^4
#     return (sig + bkg) * bin_width
# end

likelihood = let h = hist_data, centers = bin_centers
    observed_counts = h.weights

    logfuncdensity(function (p)

        function bin_ll(i)
            μ = fit_function(p, centers[i])
            if !isfinite(μ) || μ <= 0
                return -Inf
            end
            return logpdf(Poisson(μ), observed_counts[i])
        end

        ll = 0.0
        for i in eachindex(observed_counts)
            val = bin_ll(i)
            if !isfinite(val)
                return -Inf
            end
            ll += val
        end
        return ll
    end)
end

prior = distprod(
    N_sig = Uniform(0.0, 1000.0), #has to start with 0, 250.17
    p0 = Uniform(0.0, 5000.0),                        # 1852.67
    p1 = Uniform(-3000.0, 3000.0),                    # -1287.07
    p2 = Uniform(-2000.0, 2000.0),                    # 564.62
    p3 = Uniform(-1000.0, 1000.0),                    # -118.19
    p4 = Uniform(-1000.0, 1000.0)                     # -107.38
)

posterior = PosteriorMeasure(likelihood, prior)

samples = bat_sample(
    posterior,
    TransformedMCMC(
        proposal = RandomWalk(),
        nsteps = Int(1e5),
        nchains = 4
    )
).result

println("Mode: ", mode(samples))
println("Mean: ", mean(samples))
println("Std:  ", std(samples))

best = mode(samples)

N_sig_samples = getproperty.(samples.v, :N_sig)
println("N_sig mean: ", mean(N_sig_samples))
println("N_sig std: ", std(N_sig_samples))
println("N_sig 95% CI: ", quantile(N_sig_samples, [0.025, 0.975]))

yerr = sqrt.(Float64.(n))

bkg_counts = [begin
    t = (x - x_mid) / x_half
    (best.p0 + best.p1*t + best.p2*t^2 + best.p3*t^3 + best.p4*t^4) * bin_width
end for x in bin_centers]

fit_counts = [fit_function(best, x) for x in bin_centers]

fit_minus_bkg = fit_counts .- bkg_counts
data_minus_bkg = Float64.(n) .- bkg_counts

p_top = plot(bin_centers, n,
    seriestype = :scatter,
    yerr = yerr,
    label = "Data",
    ylabel = "Events",
    legend = :topright
)

plot!(p_top, bin_centers, x -> fit_function(best, x),
    label = "Fit",
    linewidth = 2
)

plot!(p_top, bin_centers,
    x -> begin
        idx = clamp(Int(floor((x - x_min) / bin_width)) + 1, 1, length(s_pdf))
        best.N_sig * s_pdf[idx] * bin_width
    end,
    label="Signal",
    linestyle=:dash,
    linewidth=2
)

plot!(p_top, bin_centers,
    x -> begin
        t = (x - x_mid) / x_half
        (best.p0 + best.p1*t + best.p2*t^2 + best.p3*t^3 + best.p4*t^4) * bin_width
    end,
    label="Background",
    linestyle=:dot,
    linewidth=2
)

plot!(p_top, xlabel="")

p_bottom = plot(bin_centers, data_minus_bkg,
    seriestype = :scatter,
    yerr = yerr,
    label = "Data - Bkg",
    xlabel = "m_γγ (GeV)",
    ylabel = "Events - Bkg",
    legend = :topright
)

plot!(p_bottom, bin_centers, fit_minus_bkg,
    label="Fit - Bkg",
    linewidth=2
)

hline!(p_bottom, [0.0],
    linestyle=:dash,
    linewidth=2,
    label="Background"
)

p_combined = plot(
    p_top,
    p_bottom,
    layout = @layout([a; b]),
    size = (800, 800),
    link = :x
)

display(p_combined)

findmode_result = bat_findmode(
    posterior,
    OptimAlg(optalg = Optim.NelderMead(), init = ExplicitInit([mode(samples)]))
)

fit_par_values = findmode_result.result
println("Refined mode: ", fit_par_values)

pC = plot(bin_centers, (p, x) -> fit_function(p, x), samples,
    xlabel = "m_γγ (GeV)",
    ylabel = "Events",
    title = "Data, Model and Best Fit",
    linewidth = 2,
    framestyle = :box,
    legend = :topright
)

plot!(pC, bin_centers, n,
    seriestype = :scatter,
    yerr = yerr,
    label = "Data",
    color = :blue
)

plot!(pC, bin_centers,
    x -> fit_function(fit_par_values, x),
    label = "Best Fit",
    color = :black,
    linestyle = :dot,
    linewidth = 2
)

display(pC)

pPost = plot(samples)
display(pPost)


#-----code not included in our original fit function:------
using Plots.PlotMeasures

# =========================
# BAYESIAN STATS — N_sig
# =========================
N_sig_samples = [s.N_sig for s in samples.v]

N_sig_mean = mean(N_sig_samples)
N_sig_std = std(N_sig_samples)
N_sig_median = median(N_sig_samples)

N_sig_lower_68 = quantile(N_sig_samples, 0.16)
N_sig_upper_68 = quantile(N_sig_samples, 0.84)
N_sig_lower_95 = quantile(N_sig_samples, 0.025)
N_sig_upper_95 = quantile(N_sig_samples, 0.975)

# =========================
# POSTERIOR PDF OF N_sig
# =========================
# pdf_posterior = histogram(N_sig_samples,
#     bins = 60,
#     label = "P(N_sig | data)",
#     xlabel = "Signal yield (N_sig)",
#     ylabel = "Count",
#     title = "Posterior: Signal Yield N_sig",
#     color = :turquoise4,
#     framestyle = :box,
#     size = (900, 600),
#     fillalpha = 1.0,
#     left_margin = 5mm,
#     bottom_margin = 3mm,
#     right_margin = 3mm,
#     top_margin = 3mm
# )

# vline!([0.0], label="No signal (N_sig=0)", lw=2, color=:red,    linestyle=:dash)
# vline!([N_sig_mean], label="Mean=$(round(N_sig_mean, digits=1))", lw=2, color=:yellow, linestyle=:dash)
# vspan!([N_sig_lower_68, N_sig_upper_68], alpha=0.2, color=:steelblue, label="68% CI")

# display(pdf_posterior)

# =========================
# BEST FIT PLOT
# =========================
bkg_counts = [begin
    t = (x - x_mid) / x_half
    (best.p0 + best.p1*t + best.p2*t^2 + best.p3*t^3 + best.p4*t^4) * bin_width
end for x in bin_centers]

sig_counts = [begin
    idx = clamp(Int(floor((x - x_min) / bin_width)) + 1, 1, length(s_pdf))
    best.N_sig * s_pdf[idx] * bin_width
end for x in bin_centers]

expected_bestfit = sig_counts .+ bkg_counts

# p_fit = plot(bin_centers, n,
#     seriestype = :scatter,
#     yerr = sqrt.(Float64.(n)),
#     label = "Data",
#     xlabel = "m_γγ (GeV)",
#     ylabel = "Events / $(bin_width) GeV",
#     title = "H → γγ: Bayesian Fit",
#     framestyle = :box,
#     grid = false,
#     legend = :topright,
#     markersize = 5,
#     color = :black,
#     size = (900, 600),
#     left_margin = 5mm,
#     bottom_margin = 3mm,
#     right_margin = 2mm
# )

# # Background
# bar!(bin_centers, bkg_counts,
#     bar_width = bin_width * 0.9,
#     label = "Background (polynomial)",
#     color = :turquoise4,
#     alpha = 0.5
# )

# # Signal stacked on top of background
# bar!(bin_centers, sig_counts,
#     bar_width = bin_width * 0.9,
#     bottom = bkg_counts,
#     label = "Signal (N_sig=$(round(best.N_sig, digits=1)))",
#     color = :red,
#     alpha = 0.7
# )

# # Total fit line
# plot!(bin_centers, expected_bestfit,
#     label = "Total fit (S+B)",
#     lw = 2,
#     color = :black,
#     linestyle = :solid,
#     marker = :circle,
#     markersize = 3
# )

# # Data on top
# scatter!(bin_centers, Float64.(n),
#     yerr = sqrt.(Float64.(n)),
#     label = "",
#     color = :black,
#     markersize = 7
# )

# display(p_fit)

# ── Styled two-panel figure ───────────────────────────────────────────────────

# 1. Main panel: BAT credible band as base (gives the colored bands + median)
p_main = plot(
    bin_centers,
    (p, x) -> fit_function(p, x),
    samples;
    xlabel = "",
    ylabel = "Events / $(bin_width) GeV",
    title = "H → γγ: Bayesian Fit",
    linewidth = 2,
    framestyle = :box,
    legend = :topright,
    legendfontsize = 12,
    grid = true,
    gridalpha = 0.3,
    minorgrid = false,
    minorgridalpha = 0.15
)

# Data points — blue filled circles matching pC style
plot!(p_main, bin_centers, n;
    seriestype = :scatter,
    yerr = yerr,
    label = "Data",
    color = :blue,
    markersize = 5,
    markerstrokecolor = :black,
    
)

# Best-fit total only — no separate signal/bkg component lines
plot!(p_main, bin_centers, x -> fit_function(best, x);
    label = "Best fit (S+B)",
    color = :black,
    linestyle = :dash,
    linewidth = 2,
)

# 2. Residuals strip
p_resid = plot(bin_centers, data_minus_bkg;
    seriestype = :scatter,
    yerr = yerr,
    label = "Data − Bkg",
    color = :black,
    markersize = 4,
    markerstrokecolor = :black,
    xlabel = "m_γγ (GeV)",
    ylabel = "Residuals",
    framestyle = :box,
    legend = :topright,
    legendfontsize = 10,
    grid = true,
    gridalpha = 0.3,
    minorgrid = false,
    minorgridalpha = 0.15
)

plot!(p_resid, bin_centers, fit_minus_bkg;
    label = "Fit − Bkg",
    color = :tomato,
    linewidth = 2,
)

hline!(p_resid, [0.0];
    linestyle = :dash,
    linewidth = 1,
    color = :gray,
    label = "",
)

# 3. Combine: 70/30 height split, shared x-axis
p_combined = plot(
    p_main,
    p_resid;
    layout = @layout([a{0.7h}; b{0.3h}]),
    size = (800, 700),
    link = :x,
)

display(p_combined)

end