
using Random, Chairmarks, CSV, DataFrames, Statistics, StatsBase

include("ebus/EBUS_bench.jl")
include("bus_julia/BUS_optimized_bench.jl")
include("alias_table/ALIAS_TABLE_bench.jl")

median_time(b) = median([x.time for x in b.samples])

rng = Xoshiro(42)

Ns = [10^i for i in 3:7]

ts_static = Dict(
  "EBUS" => Float64[],
  "BUS_jl" => Float64[],
  "ALIAS_TABLE" => Float64[]
)
for N in Ns
	w = initialize_weights_EBUS(rng, FixedSizeWeightVector, N)
	t_static_EBUS = 10^9 * median_time(@be static_samples_EBUS($rng, $w, $N) seconds=20) / N
	push!(ts_static["EBUS"], t_static_EBUS)

    ds = initialize_sampler_BUS_opt(rng, N)
    t_static_BUS_opt = 10^9 * median_time(@be static_samples_BUS_opt($rng, $ds, $N) seconds=20) / N
    push!(ts_static["BUS_jl"], t_static_BUS_opt)

    al = initialize_ALIAS_TABLE(rng, N)
    t_static_AL = 10^9 * median_time(@be static_samples_ALIAS_TABLE($rng, $al, $N) seconds=20) / N
    push!(ts_static["ALIAS_TABLE"], t_static_AL)
end
df = DataFrame(ts_static)
file = "data/static.csv"
CSV.write(file, df; append = isfile(file), writeheader = true)

ts_dynamic_fixed_dom = Dict(
  "EBUS" => Float64[],
  "BUS_jl" => Float64[]
)
for N in Ns
	w = initialize_weights_EBUS(rng, FixedSizeWeightVector, N)
	t_dynamic_fixed_dom_EBUS = 10^9 * median_time(@be dynamic_samples_fixed_dom_EBUS($rng, $w, $N) seconds=20) / N
	push!(ts_dynamic_fixed_dom["EBUS"], t_dynamic_fixed_dom_EBUS)

    ds = initialize_sampler_BUS_opt(rng, N)
    t_dynamic_fixed_dom_BUS_opt = 10^9 * median_time(@be dynamic_samples_fixed_dom_BUS_opt($rng, $ds, $N) seconds=20) / N
    push!(ts_dynamic_fixed_dom["BUS_jl"], t_dynamic_fixed_dom_BUS_opt)
end
df = DataFrame(ts_dynamic_fixed_dom)
file = "data/dynamic_fixed.csv"
CSV.write(file, df; append = isfile(file), writeheader = true)

ts_dynamic_var_dom = Dict(
  "EBUS" => Float64[],
  "BUS_jl" => Float64[]
)
for N in Ns
	insertion_order = N .+ randperm(rng, 9*N)
	steps = length(insertion_order)
	t_dynamic_var_dom_EBUS = @be initialize_weights_EBUS(rng, WeightVector, N) dynamic_samples_variable_dom_EBUS($rng, _, $insertion_order) evals=1 seconds=20
	t_dynamic_var_dom_EBUS = median_time(t_dynamic_var_dom_EBUS)
	t_dynamic_var_dom_EBUS *= 10^9 / steps
	push!(ts_dynamic_var_dom["EBUS"], t_dynamic_var_dom_EBUS)

        t_dynamic_var_dom_BUS_opt = @be initialize_sampler_BUS_opt(rng, N) dynamic_samples_variable_dom_BUS_opt($rng, _, $insertion_order) evals=1 seconds=20
    t_dynamic_var_dom_BUS_opt = median_time(t_dynamic_var_dom_BUS_opt)
        t_dynamic_var_dom_BUS_opt *= 10^9 / steps
    push!(ts_dynamic_var_dom["BUS_jl"], t_dynamic_var_dom_BUS_opt)
end
df = DataFrame(ts_dynamic_var_dom)
file = "data/dynamic_variable.csv"
CSV.write(file, df; append = isfile(file), writeheader = true)

ts_dynamic_dec_dom = Dict(
  "EBUS" => Float64[],
  "BUS_jl" => Float64[]
)
for N in Ns
    perm = randperm(rng, N)
    steps = N - fld(N, 10)

    t_dynamic_dec_dom_EBUS = @be initialize_weights_EBUS(rng, FixedSizeWeightVector, N) dynamic_samples_decreasing_dom_EBUS($rng, _, $perm) evals=1 seconds=20
	t_dynamic_dec_dom_EBUS = median_time(t_dynamic_dec_dom_EBUS)
	t_dynamic_dec_dom_EBUS *= 10^9 / steps
	push!(ts_dynamic_dec_dom["EBUS"], t_dynamic_dec_dom_EBUS)

    t_dynamic_dec_dom_BUS_opt = @be initialize_sampler_BUS_opt(rng, N) dynamic_samples_decreasing_dom_BUS_opt($rng, _, $perm) evals=1 seconds=20
    t_dynamic_dec_dom_BUS_opt = median_time(t_dynamic_dec_dom_BUS_opt)
    t_dynamic_dec_dom_BUS_opt *= 10^9 / steps
    push!(ts_dynamic_dec_dom["BUS_jl"], t_dynamic_dec_dom_BUS_opt)
end
df = DataFrame(ts_dynamic_dec_dom)
file = "data/dynamic_decreasing.csv"
CSV.write(file, df; append = isfile(file), writeheader = true)

function decaying_weights_sampling_exact(n, t)
    w = FixedSizeWeightVector(n)
    for i in 1:n
        w[i] = (2.0 + (1/(100*n)) * i)^1000
    end
    for _ in 1:t
        for i in 1:n
            w[i] /= 2.0 + (1/(100*n)) * i
        end
    end
    c = countmap(rand(w) for _ in 1:1000000)
    allvals = 1:n
    return [(x in keys(c) ? c[x] : 0) for x in allvals] ./ 1000000
end

open("data/ebus_numerical.csv", "w") do io
    for t in 1:500
        probs = decaying_weights_sampling_exact(100, t)
        println(io, join(probs, ","))
    end
end