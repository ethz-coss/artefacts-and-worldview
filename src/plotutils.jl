using StatsBase, CategoricalArrays, CairoMakie

romanenum = ["(i)", "(ii)", "(iii)", "(iv)", "(v)", "(vi)", "(vii)"]

function plot_artefacts(artefacts, position; axargs...)
    if isempty(artefacts)
        ax = Axis(position; axargs...)
        hidespines!(ax); hidedecorations!(ax)
    else
        ax = Axis(position; xticks=artefacts, xtickcolor=:red, bottomspinecolor=:red, xtrimspine=true, xtickalign=1, yticksvisible=false, axargs...)
        hidespines!(ax, :l, :r, :t); hidedecorations!(ax, ticks=false)
        xpos = (artefacts[1] + artefacts[end]) / 2
        text!(ax, xpos, .1, text=rich("artefacts", font=:italic), align = (:center, :bottom), color=:red, fontsize=10)
    end
    ax
end


function heatmap_facet(df; position, x, y, N, heatmapargs, color=:red, axkwargs...)
    ax = Axis(position; axkwargs...)
    dftv = combine(groupby(df, [y, x]), nrow)
    heatmap!(ax, dftv[:, x], dftv[:, y], dftv.nrow ./ N; heatmapargs...)
    
    dfmv = combine(groupby(df, y), x .=> [mean, iqr, median] .=> [:mean, :iqr, :median])
    scatter!(ax, dfmv[:, :mean], dfmv[:, y]; color=color, markersize=4)
    scatter!(ax, dfmv[:, :median], dfmv[:, y]; color=color, marker=:vline, markersize=10)

    grmask = dfmv[:, :iqr] .>.05

    rangebars!(ax, dfmv[grmask, y], dfmv[grmask, :mean] .- .05, dfmv[grmask, :mean] .- dfmv[grmask, :iqr]./2, color = color, direction=:x, linewidth=.5)
    rangebars!(ax, dfmv[grmask, y], dfmv[grmask, :mean] .+ .05, dfmv[grmask, :mean] .+ dfmv[grmask, :iqr]./2, color = color, direction=:x, linewidth=.5)

    rangebars!(ax, dfmv[.!grmask, y], dfmv[.!grmask, :mean] .- dfmv[.!grmask, :iqr]./2, dfmv[.!grmask, :mean] .+ dfmv[.!grmask, :iqr]./2, color = color, direction=:x, linewidth=.5)

    ax
end


function simheatmap!(
    df;
    x, y,
    row, col,
    roworder=unique(df[:, row]),
    colorder=unique(df[:, col]),
    artefacts,
    axargs = (
        xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false,
        ytickformat="{:.1f}",
        xticks=([0,.5,1], ["0", "0.5", "1"]), xminorticks=0:.1:1, xminorticksvisible=true,
        ylabel="ϵ", ylabelrotation = 0,
    ),
    N=combine(groupby(df, [row, col, y]), nrow)[1, :nrow],
    heatmapargs = (colorrange=(0,1), colormap = :YlGnBu_3),
    fig = Figure(size=(800, 400)),
    titles = colorder,
)

    alayout = GridLayout(fig[1, :]); xlayout = GridLayout(fig[2, :])
    df_grouped = groupby(df, [col, row])

    axes = Array{Any}(undef, length(roworder)+1, length(colorder))
    for (j,(c,title)) = enumerate(zip(colorder, titles))
        axes[1,j] = plot_artefacts(get(artefacts, c, Float64[]), alayout[1, j]; title)
        for (i,r) = enumerate(roworder)
            axes[i+1,j] = heatmap_facet(df_grouped[Dict(row => r, col => c)]; position=xlayout[i, j], x, y, N, heatmapargs, axargs...)
        end
    end

    linkaxes!(axes...)
    for ax in axes[2:end-1, :]; ax.xticklabelsvisible=false; end
    for ax in axes[:, 2:end]; ax.yticklabelsvisible=false; ax.ylabelvisible=false; end

    ylabellayout = GridLayout(fig[2, end+1], tellheight=true)
    Label(ylabellayout[1,:], L"$α^t \to \hat{α}$", rotation = pi/2, tellheight=false)
    Label(ylabellayout[2,:], L"$β^t \to \hat{β}$", rotation = pi/2, tellheight=false)

    Colorbar(fig.layout[2, end+1]; 
        spinewidth=0, label="% of population", labelpadding=-20, ticks=([0, .5, 1], ["0", "", "100"]), minorticks=0:.1:1, minorticksvisible=true, heatmapargs...)
    
    rowsize!(fig.layout, 1, Relative(.1))
    rowgap!(fig.layout, 1)
    colgap!(fig.layout, 1, 1)
    fig
end

function ripley_plot!(
    df;
    col, colorder=unique(df[:, col]), 
    hue, hueorder=unique(df[:, hue]),
    fig = Figure(size=(800, 300), 
    title="Ripley's L"),
    labels=hueorder,
    colors=Makie.wong_colors(),
    axargs = (
        xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false,
        ytickformat="{:.2f}", xtickformat="{:.1f}", 
        xlabel="radii",
        yticks = -.5:.25:.5, 
        limits=(extrema(radii), (-.25, .5)),
        ylabel="H", ylabelrotation = 0,
    ))
    df_grouped = groupby(df, [col, hue])
    axes = Array{Any}(undef, 1, length(colorder))
    for (j,c) = enumerate(colorder)
        axes[j] = ax = Axis(fig[:,j]; title = "ϵ=$(c)", axargs...)
        hlines!(ax, 0; color=:black, linestyle=:dash, linewidth=1)
        for (h, color, label) in zip(hueorder, colors, labels)
            lines!(ax, radii, Vector(df_grouped[Dict(hue => h, col => c)][1, Not([hue, col])]); color, label)
        end
    end
    linkaxes!(axes...)
    for ax in axes[:, 2:end]; ax.yticklabelsvisible=false; ax.ylabelvisible=false; end
    fig[:, end+1] = Legend(fig, axes[end], framevisible = false, tellheight=false; tellwidth = true)
    fig
end

function qqfacet(
    df;
    x, 
    color,
    row, col,
    roworder=unique(df[:, row]),
    colorder=unique(df[:, col]),
    axargs = (
        xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false,
        limits=(nothing, nothing, -.5, 1.5)
    ),
    fig = Figure(size=(800, 400)),
    titles = colorder,
)

    xlayout = GridLayout(fig[:, :])
    df_grouped = groupby(df, [col, row])

    axes = Array{Any}(undef, length(roworder), length(colorder))
    for (j,(c,title)) = enumerate(zip(colorder, titles))
        for (i,r) = enumerate(roworder)
            dg = df_grouped[Dict(row => r, col => c)]
            ax = axes[i,j] =Axis(xlayout[i,j]; title, axargs...)
            for cc in color
                qqnorm!(ax, dg[dg[:, color] .== :cc, x], qqline = :fitrobust, markersize=4, label=cc)
            end
        end
    end
    linkaxes!(axes...)
    for ax in axes[1:end-1, :]; ax.xticklabelsvisible=false; end
    for ax in axes[:, 2:end]; ax.yticklabelsvisible=false; ax.ylabelvisible=false; end

    fig
end