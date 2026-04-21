
using WeightVectors

function initialize_weights_EBUS(rng, WType, N)
    w = WType(N)
    for i in eachindex(w)
	w[i] = abs(randn(rng))
    end
    return w
end

static_samples_EBUS(rng, ds, N) = rand(rng, ds, N)

function dynamic_samples_fixed_dom_EBUS(rng, w, N)
    s = Vector{Int}(undef, N)
    @inbounds for i in 1:N
	s[i] = rand(rng, w)
        w[rand(rng, 1:N)] = abs(randn(rng))
    end
    return s
end

function dynamic_samples_variable_dom_EBUS(rng, w, insertion_order)
    steps = length(insertion_order)
    s = Vector{Int}(undef, steps)
    resize!(w, length(w) + steps)
    @inbounds for i in 1:steps
        s[i] = rand(rng, w)
        w[insertion_order[i]] = abs(randn(rng))
    end
    return s
end

function dynamic_samples_decreasing_dom_EBUS(rng, w, removal_order)
    initial_n = length(removal_order)
    steps = initial_n - fld(initial_n, 10)
    s = Vector{Int}(undef, steps)
    @inbounds for i in 1:steps
        s[i] = rand(rng, w)
        w[removal_order[i]] = 0.0
    end
    return s
end
