module ToyGeneration

using Distributions, Plots, Measures
using ..HiggsAnalysis3

export toys

# ===== Toy Experiments: Generate signal + background toys =====

function generate_toy_experiment(signal, background, mu)
    # Generate Poisson counts for each bin based on signal + background
    expected_counts = mu .* signal .+ background
    toy_counts = [rand(Poisson(λ)) for λ in expected_counts]
    return toy_counts
end

# Generate 100 toy experiments using fit window data
toys = [generate_toy_experiment(fit_signal, fit_background, best_fit.mu) for _ in 1:100]

# Plot toy experiments
p_toys = plot(
    xlabel = "m4l [GeV]",
    ylabel = "Events / 2 GeV",
    title = "Toy Experiments: Signal + Background",
    framestyle = :box,
    size = (900, 600),
    legend = :topright,
    left_margin = 5mm,
    bottom_margin = 3mm,
    right_margin = 2mm
)

# Show all 10 toy experiments
for (i, toy) in enumerate(toys)
    bar!(p_toys, fit_bin_centers, toy, 
          label = "Toy $i", 
          alpha = 0.6,)
end

# Overlay the expected signal + background (μ=1)
expected_nominal = fit_signal .+ fit_background
plot!(p_toys, fit_bin_centers, expected_nominal,
      label = "Expected (μ=1, S+B)",
      seriestype = :steppost,
      color = :green,
      lw = 3,
      linestyle = :dash)

# Overlay best fit expectation
expected_bestfit = best_fit.mu .* fit_signal .+ fit_background
plot!(p_toys, fit_bin_centers, expected_bestfit,
      label = "Best Fit (μ=$(round(best_fit.mu, digits=2)), μS+B)",
      seriestype = :steppost,
      color = :red,
      lw = 3)

display(p_toys)

end