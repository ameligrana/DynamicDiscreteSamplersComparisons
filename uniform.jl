
using Random
using BenchmarkTools
using Plots

function benchmark_uniform_sampling(; seconds=2.0)
    sizes = 10 .^ (3:7)

    medians_ns = Float64[]

    for n in sizes
        v = collect(1:n)
        GC.gc()
        trial = @benchmark rand($v) seconds=seconds
        est = median(trial)
        push!(medians_ns, est.time)
    end

    p = plot(
        sizes,
        medians_ns;
        marker = :circle,
        linewidth = 2,
        xlabel = "Vector length",
        ylabel = "Median sampling time (ns)",
        title = "Uniform random sampling from Vector",
        xscale = :log10,
        xticks = sizes,
        legend = false,
        grid = true,
    )

    return p
end

p = benchmark_uniform_sampling(seconds=2.0)
display(p)
savefig(p, "figures/uniform_sampling_benchmark.png")