
using DynamicSampling

function initialize_sampler_BUS_opt(rng, N)
	ds = DynamicSampler(rng)
	append!(ds, 1:N, (abs(randn(rng)) for _ in 1:N))
	return ds
end

static_samples_BUS_opt(rng, ds, N) = rand(ds, N)

function dynamic_samples_fixed_dom_BUS_opt(rng, ds, N)
	s = Vector{Int}(undef, N)
	@inbounds for i in 1:N
		s[i] = rand(ds)
		r = rand(rng, 1:N)
		delete!(ds, r)
		push!(ds, r, abs(randn(rng)))
	end
	return s
end

function dynamic_samples_variable_dom_BUS_opt(rng, ds, insertion_order)
	steps = length(insertion_order)
	s = Vector{Int}(undef, steps)
	@inbounds for i in 1:steps
		s[i] = rand(ds)
		push!(ds, insertion_order[i], abs(randn(rng)))
	end
	return s
end

function dynamic_samples_decreasing_dom_BUS_opt(rng, ds, removal_order)
	initial_n = length(removal_order)
	steps = initial_n - fld(initial_n, 10)
	s = Vector{Int}(undef, steps)
	@inbounds for i in 1:steps
		s[i] = rand(ds)
		delete!(ds, removal_order[i])
	end
	return s
end

