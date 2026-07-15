module UnweightedPowerFit

using DelimitedFiles, Plots, LaTeXStrings, LsqFit

export fit_unweighted_power

function fit_unweighted_power(; csv_path="sensitivity_results.csv")
    # 1. Load Data
    if !isfile(csv_path)
        error("File '$csv_path' not found.")
    end
    data = readdlm(csv_path, ',', Float64, skipstart=1)
    
    mu      = data[:, 1]
    medians = data[:, 2]
    q16     = data[:, 3]
    q84     = data[:, 4]

    # 2. Define Model: y = C * μ^α + K
    model(x, p) = p[1] .* x.^p[2] .+ p[3]
    
    # Initial Guess [C, alpha, K]
    p0 = [2.5, 2.0, 0.0]

    # 3. Perform UNWEIGHTED fit
    # By omitting the 'w' argument, curve_fit treats all points equally
    fit = curve_fit(model, mu, medians, p0)
    C, alpha, K = fit.param

    println("--- Unweighted Power Law Results ---")
    println("Formula: y = $(round(C, digits=3)) * μ^$(round(alpha, digits=3)) + $(round(K, digits=3))")

    # 4. Plotting
    p = scatter(mu, medians, yerror=(medians.-q16, q84.-medians),
                label="Toy Medians", color=:blue, markersize=5,
                xlabel=L"\mu", ylabel=L"\mathrm{median}(\log BF)",
                title="Median Log Bayes Factor vs. Signal Strength for 4l Channel",
                framestyle=:box, size=(800, 500))

    mu_fine = range(0.01, stop=maximum(mu), length=200)
    plot!(mu_fine, model(mu_fine, fit.param), lw=3, color=:darkred,
          label=L"Fit: %$(round(C,digits=2)) \mu^{%$(round(alpha,digits=2))} + %$(round(K,digits=2))")

    hline!([4.5],  color=:orange, ls=:dash, lw=1.5, label=L"3\sigma\ (\log BF=4.5)")
    hline!([12.5], color=:green,  ls=:dash, lw=1.5, label=L"5\sigma\ (\log BF=12.5)")

    display(p)
    return fit.param
end

fit_unweighted_power(; csv_path="sensitivity_results.csv")

end # module