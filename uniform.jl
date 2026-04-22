
using Random
using BenchmarkTools
using Plots

function benchmark_uniform_sampling()
    sizes = 10 .^ (3:8)

    medians_ns = Float64[]
    rng = Xoshiro(42)

    vecs = Vector{Int}[]
    for n in sizes
        v = collect(1:n)
        push!(vecs, v)
    end
    GC.gc()
    for v in vecs
        trial = @benchmark rand($rng, $v)
        est = median(trial)
        push!(medians_ns, est.time)
    end

    p = plot(sizes, medians_ns;
        marker = :circle, linewidth = 2,
        xlabel = "Vector length", ylabel = "Median sampling time (ns)",
        title = "Uniform random sampling from Vector", xscale = :log10,
        xticks = sizes, legend = false, grid = true)

    return p
end

p = benchmark_uniform_sampling()
display(p)
savefig(p, "figures/uniform_sampling_benchmark.png")