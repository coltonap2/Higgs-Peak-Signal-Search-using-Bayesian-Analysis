module diphoton_toy_generation

using Distributions
using Plots
using Measures

using ..diphoton_fit

export toys, generate_toy_experiment  # can't leave export blank

function generate_toy_experiment(best, s_pdf, bin_centers, bin_width, x_min, x_max)
    x_mid  = (x_min + x_max) / 2
    x_half = (x_max - x_min) / 2

    expected_counts = map(enumerate(bin_centers)) do (i, x)
        t   = (x - x_mid) / x_half
        sig = best.N_sig * s_pdf[i]
        bkg = best.p0 + best.p1*t + best.p2*t^2 + best.p3*t^3 + best.p4*t^4
        return (sig + bkg) * bin_width
    end

    toy_counts = [rand(Poisson(lamby_bamby)) for lamby_bamby in expected_counts]
    return toy_counts
end

# best comes from diphoton_fit (your mode/findmode result — make sure it's exported)
toys = [generate_toy_experiment(best, s_pdf, bin_centers, bin_width, x_min, x_max) for _ in 1:300]

# Plot
p_toys = plot(
    xlabel = "m_γγ [GeV]",      # fixed: was m4l which is the golden channel label
    ylabel = "Events / 2 GeV",
    title = "Diphoton Toy Experiments: Signal + Background",
    framestyle = :box,
    size = (900, 600),
    legend = :topright,
    left_margin = 5mm,
    bottom_margin = 3mm,
    right_margin = 2mm
)

for (i, toy) in enumerate(toys)
    bar!(p_toys, bin_centers, toy,   # fixed: was fit_bin_centers, which doesn't exist here
          label = "Toy $i", 
          alpha = 0.6)
end

# Overlay the nominal expected counts (using best-fit parameters as-is)
# In diphoton there is no separate fit_signal/fit_background vector to pull from,
# so we compute expected counts the same way the generator does
expected_nominal = map(enumerate(bin_centers)) do (i, x)
    x_mid  = (x_min + x_max) / 2
    x_half = (x_max - x_min) / 2
    t   = (x - x_mid) / x_half
    sig = best.N_sig * s_pdf[i]
    bkg = best.p0 + best.p1*t + best.p2*t^2 + best.p3*t^3 + best.p4*t^4
    return (sig + bkg) * bin_width
end

plot!(p_toys, bin_centers, expected_nominal,
      label = "Expected (best fit S+B)",
      seriestype = :steppost,
      color = :green,
      lw = 3,
      linestyle = :dash)

display(p_toys)

end