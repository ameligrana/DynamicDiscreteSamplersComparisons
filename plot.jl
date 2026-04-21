
using Plots, CSV, DataFrames

p = nothing

function subset_and_rename(df, cols, rename_map)
    df_cols = propertynames(df)
    present_cols = filter(c -> c in df_cols, cols)
    df = select(df, present_cols)
    present_rename_map = Dict(k => v for (k, v) in rename_map if k in present_cols)
    rename!(df, present_rename_map)
    return df
end

ts_static = CSV.read("data/static.csv", DataFrame)
ts_dynamic_fixed_dom = CSV.read("data/dynamic_fixed.csv", DataFrame)
ts_dynamic_var_dom = CSV.read("data/dynamic_variable.csv", DataFrame)
ts_dynamic_dec_dom = CSV.read("data/dynamic_decreasing.csv", DataFrame)

rename_map = Dict(:EBUS => :EBUS, :FOREST_OF_TREES => :FT, :PROPOSAL_ARRAY => :DPA, :BUS => :BUS)
ts_static = subset_and_rename(ts_static, [:EBUS, :BUS, :FOREST_OF_TREES, :PROPOSAL_ARRAY], rename_map)

ts_dynamic_fixed_dom = subset_and_rename(ts_dynamic_fixed_dom, [:EBUS, :BUS, :FOREST_OF_TREES, :PROPOSAL_ARRAY], rename_map)

rename_map = Dict(:EBUS => :EBUS, :PROPOSAL_ARRAY => :DPA, :BUS => :BUS)
ts_dynamic_var_dom = subset_and_rename(ts_dynamic_var_dom, [:EBUS, :BUS, :PROPOSAL_ARRAY], rename_map)

rename_map = Dict(:EBUS => :EBUS, :FOREST_OF_TREES => :FT, :PROPOSAL_ARRAY => :DPA, :BUS => :BUS)
ts_dynamic_dec_dom = subset_and_rename(ts_dynamic_dec_dom, [:EBUS, :BUS, :FOREST_OF_TREES, :PROPOSAL_ARRAY], rename_map)

markers = Dict(:EBUS => :circle, :FT => :rect, :DPA => :utriangle, :BUS => :xcross)
colors = Dict(:EBUS => 4, :FT => 2, :DPA => 1, :BUS => 3)

Ns = [10^i for i in 3:7]
perf_yticks = 0:200:800
perf_ylims = (0, 800)
subplot_title_fontsize = 14
axis_guide_fontsize = 11
axis_tick_fontsize = 10
legend_fontsize = 11
figure_title_fontsize = 17

single_title_fontsize = 10
single_axis_guide_fontsize = 8
single_axis_tick_fontsize = 7
single_legend_fontsize = 7

function plot_vals(Ns, ts, order, title, xlabel, ylabel, pname; show_legend=true, save_outputs=true,
    show_xlabel=true, show_ylabel=true)
    p = nothing
    present_order = filter(k -> k in propertynames(ts), order)
    for (i, k) in enumerate(present_order)
        if i == 1
            p = plot(Ns, ts[!, k], xscale=:log10, marker=markers[k], xticks=Ns, 
                 xlabel=(show_xlabel ? xlabel : ""), ylabel=(show_ylabel ? ylabel : ""), markersize=6, line = (2, :dash),
                 title=title, label=(k == :DPA ? (string(k) * "*") : string(k)),
                ylims=perf_ylims, yticks=perf_yticks, widen = true, legend=(show_legend ? :topleft : false),
                                    titlefontsize=subplot_title_fontsize, guidefontsize=axis_guide_fontsize,
                                    tickfontsize=axis_tick_fontsize, legendfontsize=legend_fontsize,
                                    right_margin=10Plots.mm, color = colors[k], dpi=1000)
        else
            plot!(Ns, ts[!, k], marker=markers[k], label=(k == :DPA ? (string(k) * "*") : string(k)),
                                ylims=perf_ylims, yticks=perf_yticks, widen = true, color = colors[k], markersize=6,
                                line = (2, :dash), legendfontsize=legend_fontsize, dpi=1000)
        end
    end
    if save_outputs
        savefig(p, "figures/" * pname * ".pdf")
        savefig(p, "figures/" * pname * ".png")
    end
    return p
end

ps_static = plot_vals(Ns, ts_static, [:FT, :DPA, :BUS, :EBUS], "Static Sampling", "sampler size",
    "time per single draw (ns)", "static")

ps_dynamic_fixed = plot_vals(Ns, ts_dynamic_fixed_dom, [:FT, :DPA, :BUS, :EBUS], "Fixed Range Sampling", "sampler size",
    "time per single update & draw (ns)", "dynamic_fixed")

ps_dynamic_var = plot_vals(Ns, ts_dynamic_var_dom, [:DPA, :BUS, :EBUS], "Increasing Range Sampling", "starting sampler size",
    "time per single update & draw (ns)", "dynamic_variable")

ps_dynamic_dec = plot_vals(Ns, ts_dynamic_dec_dom, [:FT, :DPA, :BUS, :EBUS], "Decreasing Range Sampling", "starting sampler size",
    "time per single update & draw (ns)", "dynamic_decreasing")

single_plot_size = Plots.default(:size)
ps_dynamic_fixed_sub = plot_vals(Ns, ts_dynamic_fixed_dom, [:FT, :DPA, :BUS, :EBUS], "Fixed Range Sampling", "sampler size",
    "time (ns)", "dynamic_fixed"; show_legend=false, save_outputs=false,
    show_xlabel=false, show_ylabel=false)
ps_dynamic_dec_sub = plot_vals(Ns, ts_dynamic_dec_dom, [:FT, :DPA, :BUS, :EBUS], "Decreasing Range Sampling", "sampler size",
    "time (ns)", "dynamic_decreasing"; show_legend=false, save_outputs=false,
    show_xlabel=true, show_ylabel=true)
ps_dynamic_var_sub = plot_vals(Ns, ts_dynamic_var_dom, [:DPA, :BUS, :EBUS], "Increasing Range Sampling", "sampler size",
    "time (ns)", "dynamic_variable"; show_legend=false, save_outputs=false,
    show_xlabel=true, show_ylabel=false)

ps_static_sub = plot_vals(Ns, ts_static, [:FT, :DPA, :BUS, :EBUS], "Static Sampling", "sampler size",
    "time (ns)", "static"; show_legend=true, save_outputs=false,
    show_xlabel=false, show_ylabel=true)

p_perf = plot(ps_static_sub, ps_dynamic_fixed_sub, ps_dynamic_dec_sub, ps_dynamic_var_sub,
    layout=(2, 2), size=(round(Int, 2.1 * single_plot_size[1]), round(Int, 2.1 * single_plot_size[2])),
    margin=4Plots.mm,
    plot_titlefontsize=figure_title_fontsize,
    dpi=1000)
savefig(p_perf, "figures/performance_benchmarks.pdf")
savefig(p_perf, "figures/performance_benchmarks.png")

### numerical

function decaying_weights_sampling_probability(n, t)
    p = BigFloat[]
    for i in 1:n
        push!(p, BigFloat((2.0 + (1/(100*n)) * i))^(1000-t))
    end
    return Float64.(p ./ sum(p))
end

function js_divergence(a, b)
    js = Float64[]
    for (ai, bi) in zip(a, b)
    u = (ai + bi) / 2
    ta = iszero(ai) ? 0.0 : ai * log2(ai) / 2
    tb = iszero(bi) ? 0.0 : bi * log2(bi) / 2
    tu = iszero(u) ? 0.0 : u * log2(u)
    push!(js, ta + tb - tu)
   end
   return sum(js)
end

df_EBUS = CSV.read("data/ebus_numerical.csv", DataFrame, header=false)
M_EBUS = Matrix(df_EBUS)

df_FT = CSV.read("data/forest_of_trees_numerical.csv", DataFrame, header=false)
M_FT = Matrix(df_FT)

df_DPA = CSV.read("data/proposal_array_numerical.csv", DataFrame, header=false)
M_DPA = Matrix(df_DPA)

df_BUS = CSV.read("data/bus_numerical.csv", DataFrame, header=false)
M_BUS = Matrix(df_BUS)

js_EBUS = []
js_FT = []
js_DPA = []
js_BUS = []

for i in 1:100
    probs = decaying_weights_sampling_probability(100, i)
    push!(js_EBUS, js_divergence(M_EBUS[i, :], probs))
    i <= 51 && push!(js_FT, js_divergence(M_FT[i, :], probs))
    i <= 51 && push!(js_DPA, js_divergence(M_DPA[i, :], probs))
    i <= 50 && push!(js_BUS, js_divergence(M_BUS[i, :], probs))
end

p = plot(js_FT, line = (1, :dot), label="FT", 
    ylabel="divergence", xlabel="decay step",
    size=(single_plot_size[1], round(Int, 0.75 * single_plot_size[2])),
    right_margin=10Plots.mm, bottom_margin=3Plots.mm, dpi=1000,
    titlefontsize=single_title_fontsize, guidefontsize=single_axis_guide_fontsize,
    tickfontsize=single_axis_tick_fontsize, legendfontsize=single_legend_fontsize,
    legend=:topleft, color = colors[:FT])
plot!(js_DPA, line = (1, :dash), label="DPA*", color = colors[:DPA])
plot!(js_BUS, line = (1, :dashdot), label="BUS", color = colors[:BUS])
plot!(js_EBUS, line = (1, :solid), label="EBUS", color = colors[:EBUS])
savefig(p, "figures/numerical" * ".pdf")
savefig(p, "figures/numerical" * ".png")

k = 1:100
l = @layout [a b; d e]

yticks = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
p_theory = decaying_weights_sampling_probability(100, 50)
p2 = bar(k, p_theory, bar_width = 0.9, title = "EBUS", label = "theoretical",
    ylabel="probability", yticks=yticks, ylims=(0, 0.06), linecolor=2,
    color=2, fillalpha=0.25, alpha=0.25, legend=:topleft,
    titlefontsize=subplot_title_fontsize, guidefontsize=axis_guide_fontsize,
    tickfontsize=axis_tick_fontsize, legendfontsize=legend_fontsize)
bar!(p2, k, M_EBUS[50, :], fillalpha=0.25, alpha=0.25, bar_width = 0.6, label = "empirical", linecolor=1, color=1)

p3 = bar(k, p_theory, bar_width = 0.9, title = "BUS", label = "theoretical",
    yticks=yticks, ylims=(0, 0.06), linecolor=2,
    color=2, fillalpha=0.25, alpha=0.25, legend = false,
    titlefontsize=subplot_title_fontsize, guidefontsize=axis_guide_fontsize,
    tickfontsize=axis_tick_fontsize, legendfontsize=legend_fontsize)
bar!(p3, k, M_BUS[50, :], fillalpha=0.25, alpha=0.25, bar_width = 0.6, label = "empirical", linecolor=1, color=1)

p4 = bar(k, p_theory, bar_width = 0.9, title = "FT", label = "theoretical",
    xlabel="index", ylabel="probability", yticks=yticks, ylims=(0, 0.06), linecolor=2,
    color=2, fillalpha=0.25, alpha=0.25, legend = false,
    titlefontsize=subplot_title_fontsize, guidefontsize=axis_guide_fontsize,
    tickfontsize=axis_tick_fontsize, legendfontsize=legend_fontsize)
bar!(p4, k, M_FT[50, :], fillalpha=0.25, alpha=0.25, bar_width = 0.6, label = "empirical", linecolor=1, color=1)

p5 = bar(k, p_theory, bar_width = 0.9, title = "DPA*", label = "theoretical",
    xlabel="index", yticks=yticks, ylims=(0, 0.06), linecolor=2,
    color=2, fillalpha=0.25, alpha=0.25, legend = false,
    titlefontsize=subplot_title_fontsize, guidefontsize=axis_guide_fontsize,
    tickfontsize=axis_tick_fontsize, legendfontsize=legend_fontsize)
bar!(p5, k, M_DPA[50, :], fillalpha=0.25, alpha=0.25, bar_width = 0.6, label = "empirical", linecolor=1, color=1)

p = plot(p2, p3, p4, p5, layout = l,
    size=(round(Int, 2.1 * single_plot_size[1]), round(Int, 2.1 * single_plot_size[2])),
    margin=4Plots.mm,
    plot_titlevspan=0.04, dpi=1000, plot_titlefontsize=figure_title_fontsize)

savefig(p, "figures/numerical50" * ".pdf")
savefig(p, "figures/numerical50" * ".png")
