using CSV
using DataFrames
using CairoMakie
using Chain

df = DataFrame()
for f in readdir("data/outputs/"; join= true)
    bnf = basename(f)
    contains(bnf, "gsea_full") || continue
    newdf = CSV.read(f, DataFrame; stringtype=String)
    kind = contains(bnf, "_noage_") ? "noage" :
            contains(bnf, "_pred_") ? "mbage" : "gtage"
    newdf.model .= kind
    append!(df, newdf)
end

gdf = groupby(df, "model")

#-

fig = Figure(; resolution=(1200, 800));
ax1 = Axis(fig[1,1]; xlabel="model", ylabel="n significant genesets", xticks = ([1,2,3], [k.model for k in keys(gdf)]))
ax2 = Axis(fig[1,2]; xlabel="rank", ylabel="corrected p-value")
ax3 = Axis(fig[2,1:3]; xlabel="rank", ylabel="corrected p-value")

for (i, grp) in enumerate(gdf)
    barplot!(ax1, [i], [count(<(0.2), grp.qvalue)])
    grp = sort(grp, "qvalue")
    lines!(ax2, 1:nrow(grp), grp.qvalue)
    lines!(ax3, 1:nrow(grp), grp.qvalue)

end
hlines!(ax2, [0.2]; color=:black, linestyle=:dot)
hlines!(ax3, [0.2]; color=:black, linestyle=:dot)
ylims!(ax3, (0., 0.4))
xlims!(ax3, (0., 50))
Legend(fig[1,3], [MarkerElement(; color=Makie.wong_colors()[i], marker=:rect) for i in 1:3],
    [k.model for k in keys(gdf)])

fig

