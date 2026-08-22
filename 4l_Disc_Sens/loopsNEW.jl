#using Pkg
#Pkg.add(["LaTeXStrings"])
module MedianSensitivityScan

using Statistics, Random, Plots, Measures, LinearAlgebra, LaTeXStrings
using ..HiggsAnalysis3
using ..ToyGeneration
using ..ToyBATFit

export run_median_sensitivity

function run_median_sensitivity()
    # 1. Configuration
    mu_values = collect(0.0:0.1:2.5) 
    toys_per_mu = 200
    
    median_lbfs = Float64[]
    q16_lbfs = Float64[]
    q84_lbfs = Float64[]

    println("Starting Median-Based Sensitivity Scan...")

    for mu_test in mu_values
        current_mu_lbfs = Float64[]
        
        # Generate and process toys 
        for i in 1:toys_per_mu
            toy_data = generate_toy_experiment(fit_signal, fit_background, mu_test)
            lbf = get_log_bf(toy_data) 
            
            if isfinite(lbf)
                push!(current_mu_lbfs, lbf)
            end
        end

        # 2. Extract Median and 68% Quantiles 
        if !isempty(current_mu_lbfs)
            push!(median_lbfs, median(current_mu_lbfs))
            push!(q16_lbfs, quantile(current_mu_lbfs, 0.16))
            push!(q84_lbfs, quantile(current_mu_lbfs, 0.84))
        else
            push!(median_lbfs, 0.0); push!(q16_lbfs, 0.0); push!(q84_lbfs, 0.0)
        end
    end

    # 3. Quadratic Fit: f(μ) = a * μ² 
    mask = mu_values .> 0
    x_fit = mu_values[mask]
    y_fit = median_lbfs[mask]
    
    # Weights based on the 68% spread 
    sigma_fit = 0.5 .* (q84_lbfs[mask] .- q16_lbfs[mask])
    w = 1.0 ./ (max.(sigma_fit, 1e-3).^2)

    # Solve for 'a' using weighted least squares 
    a_coeff = sum(w .* (x_fit.^2) .* y_fit) / sum(w .* (x_fit.^4))
    fit_func(x) = a_coeff * x^2

    # 4. Discovery Threshold Calculation 
    # Finding μ needed for 3-sigma (log BF ≈ 4.5) and 5-sigma (log BF ≈ 12.5)
    mu_3sig = sqrt(4.5 / a_coeff)
    mu_5sig = sqrt(12.5 / a_coeff)

    println("Fitted coefficient a: ", round(a_coeff, digits=5))
    println("Expected μ for 3σ discovery: ", round(mu_3sig, digits=2))

    # 5. Visualization
    lower_err = median_lbfs .- q16_lbfs
    upper_err = q84_lbfs .- median_lbfs

    p = scatter(
        mu_values,
        median_lbfs,
        yerror = (lower_err, upper_err),
        label = L"\mathrm{Toy\ Median}",
        xlabel = L"\mathrm{True\ Signal\ Strength\ (\mu)}",
        ylabel = L"\mathrm{median}(\log BF)",
        title = "Discovery Sensitivity (Median-Based)",
        color = :blue,
        markersize = 5,
        framestyle = :box,
        size = (900, 600)
    )

    # Plot the fit line
    mu_range = range(0, stop=maximum(mu_values), length=300)
    plot!(mu_range, fit_func.(mu_range), lw=3, color=:black, label=L"a \mu^2")

    # Threshold lines
    hline!([4.5], color=:orange, linestyle=:dash, label=L"3\sigma\ (\log BF=4.5)")
    hline!([12.5], color=:green, linestyle=:dash, label=L"5\sigma\ (\log BF=12.5)")

    display(p)
    return a_coeff
end

run_median_sensitivity()

end # module