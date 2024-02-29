# Requires manuscript/Figures.jl to be run first
include("../manuscript/Figures.jl")

#-

figure2 = Figure(; size = (1100,600))

grid_fsea_dots = GridLayout(figure2[1,1])
# grid_fsea_heatmaps = GridLayout(figure2[1,2])
# grid_bug_heatmaps = GridLayout(figure2[2, 1])

gs_interval = 6
tp_interval = 1.5
dotsxlim = 3
fsea_marker_size=6

grid_fsea_latency=GridLayout(grid_fsea_dots[1,1])
ax_lat = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_fsea_latency[1,i];
				 xlabel = "z",
				 title = lab,
				 yticks = (tickrange, reverse(gssig)),
				 )
		r = Axis(grid_fsea_latency[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						  repeat(string.(3:-1:1); outer=length(tickrange)))
				 )
		(l, r)
	end
end

grid_fsea_amplitude=GridLayout(grid_fsea_dots[1,2])
ax_amp = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_fsea_amplitude[1,i];
				 xlabel = "z",
				 title = lab,
				 yticks = (tickrange, reverse(gssig)),
				 )
		r = Axis(grid_fsea_amplitude[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						  repeat(string.(3:-1:1); outer=length(tickrange)))
				 )
		(l, r)
	end
end

for axs in (ax_lat, ax_amp)
	for (i, (axl, axr)) in enumerate(axs)
		ylims!.((axl,axr), gs_interval - 2*tp_interval, gs_interval * length(gssig) + 2*tp_interval)
		xlims!(axl, -dotsxlim, dotsxlim)
		hidexdecorations!(axl; ticks=false, ticklabels=false, label=false)
		hideydecorations!(axl, ticks=i != 1, ticklabels=i != 1)
		hideydecorations!(axr, ticks=i != 3, ticklabels=i != 3)
		hidexdecorations!(axr)

		for gs in gssig
			mid = gsidx[gs] * gs_interval
			isodd(gsidx[gs]) || continue
			span = gs_interval / 2
			poly!(axl, Point2f.([(-dotsxlim, mid - span),
								 (dotsxlim, mid - span),
								 (dotsxlim, mid + span),
								 (-dotsxlim, mid + span)]);
				  color=("gray80", 0.4))
		end

	end
end
hideydecorations!(ax_amp[1][1], ticks=false)

Legend(grid_fsea_dots[2,1:2], 
	   [[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in 1:3],
		[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in 5:7]],
	   [["q < 0.01", "q < 0.1", "q < 0.2"], 
		["q < 0.2", "q < 0.1", "q < 0.01"]],
	   ["(-)", "(+)"];
	   orientation=:horizontal, tellheight=true, tellwidth=false)

#-


for (i, tp) in enumerate(tps)
    @info "$tp"
    for (j, feat) in enumerate(filter(contains("latency"), eeg_features))
        gdf = groupby(fsea_df, "eeg_feature")
        featstr = replace(feat, "peak_latency_"=>"", "_corrected"=>"c")
        @warn "$featstr"

		df =  alllms[(; eeg_feature=feat, timepoint=tp)]
		allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        for row in eachrow(subdf)
			yidx = df[!, row.geneset]
			xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
			ys = df.z[yidx]
			xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? colors_sig[cidx] : colors_sig[8 - cidx]
            violin!(ax_lat[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_lat[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
			lines!(ax_lat[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end

    for (j, feat) in enumerate(filter(contains("amp"), eeg_features))
        gdf = groupby(fsea_df, "eeg_feature")
        featstr = replace(feat, "peak_amp"=>"", "_corrected"=>"c")
        @warn "$featstr"
        df =  alllms[(; eeg_feature=feat, timepoint=tp)]
        allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        # lines!(ax_lat, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        for row in eachrow(subdf)
            # haskey(gsidx, row.geneset) || continue
            yidx = df[!, row.geneset]
            xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
            ys = df.z[yidx]
            xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
            ymed = median(ys)

            cidx = row.q₀ > 0.2  ? 4 :       # not significant
            row.q₀ < 0.01 ? 1 :       # quite significant
            row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? colors_sig[cidx] : colors_sig[8 - cidx]
            violin!(ax_amp[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_amp[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
            lines!(ax_amp[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end
    save("/home/kevin/Downloads/figure2_$tp.png", figure2)
end



############
# Figure 3 #
############

figure3 = Figure(; size = (1050, 750))


grid_future_violins = GridLayout(figure3[1,1])
grid_futfsea_dots = GridLayout(figure3[1,2])
ax_future_violins = map(enumerate((future_3m6m, future_3m12m, future_6m12m))) do (i, comm)
	ax = Axis(grid_future_violins[i,1]; ylabel = "age (months)", xticks = ([1,2], ["stool", "eeg"]),
		 xlabel=i == 3 ? "collection type" : "", title = "visit $(i == 3 ? "2" : "1") → visit $(i == 1 ? "2" : "3")")
	hidedecorations!(ax; label=false, ticks=false, ticklabels=false)
	ax
end

grid_futfsea_latency=GridLayout(grid_futfsea_dots[1,1])
ax_futlat = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_futfsea_latency[1,i];
				 xlabel = "z",
				 xticklabelsize= 10,
				 title = lab,
				 yticks = (tickrange, reverse(gssig)),
				 yticklabelsize=10,
				 )
		r = Axis(grid_futfsea_latency[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						   repeat(["v2→v3", "v1→v3", "v1→v2"]; outer=length(tickrange)))
				 )
		(l, r)
	end
end

grid_futfsea_amplitude=GridLayout(grid_futfsea_dots[1,2])
ax_futamp = let tickrange = gs_interval:gs_interval:length(gssig)*gs_interval
	map(enumerate(["N1", "P1c", "N2c"])) do (i, lab)
		l = Axis(grid_futfsea_amplitude[1,i];
				 xlabel = "z",
				 xticklabelsize=10,
				 title = lab,
				 yticks = (tickrange, reverse(gssig)),
				 yticklabelsize=10,
				 )
		r = Axis(grid_futfsea_amplitude[1,i];
				 yaxisposition = :right,
				 yticklabelsize = 6,
				 yticks = (sort([collect(tickrange);
								 collect(tickrange) .- tp_interval;
								 collect(tickrange) .+ tp_interval]),
						  repeat(["v2→v3", "v1→v3", "v1→v2"]; outer=length(tickrange)))
				 )
		(l, r)
	end
end

for axs in (ax_futlat, ax_futamp)
	for (i, (axl, axr)) in enumerate(axs)
		ylims!.((axl,axr), gs_interval - 2*tp_interval, gs_interval * length(gssig) + 2*tp_interval)
		xlims!(axl, -dotsxlim, dotsxlim)
		hidexdecorations!(axl; ticks=false, ticklabels=false, label=false)
		hideydecorations!(axl, ticks=i != 1, ticklabels=i != 1)
		hideydecorations!(axr, ticks=i != 3, ticklabels=i != 3)
		hidexdecorations!(axr)

		for gs in gssig
			mid = gsidx[gs] * gs_interval
			isodd(gsidx[gs]) || continue
			span = gs_interval / 2
			poly!(axl, Point2f.([(-dotsxlim, mid - span),
								 (dotsxlim, mid - span),
								 (dotsxlim, mid + span),
								 (-dotsxlim, mid + span)]);
				  color=("gray80", 0.4))
		end

	end
end
hideydecorations!(ax_futamp[1][1], ticks=false)

Legend(grid_future_violins[4,1], 
	   [[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in 1:3],
		[MarkerElement(; marker=:rect, color=colors_sig[i]) for i in reverse(5:7)]],
	   [["q < 0.01", "q < 0.1", "q < 0.2"], 
		reverse(["q < 0.2", "q < 0.1", "q < 0.01"])],
	   ["(-)", "(+)"];
	   orientation=:horizontal, tellheight=true, tellwidth=false, nbanks=3)

#-

for (i, comm) in enumerate((future_3m6m, future_3m12m, future_6m12m))
	ax = ax_future_violins[i]
	clrs = [colors_timepoints[i == 3 ? 2 : 1][2], colors_timepoints[i == 1 ? 2 : 3][2]]
		    
	df = DataFrame(get(comm))
	violin!(ax, repeat([1,2]; inner=size(df,1)), [df.age; df.eeg_age];
			color=repeat([(c, 0.3) for c in clrs]; inner=size(df, 1)))
	for row in eachrow(df)
		xs = [1,2] .+ rand(Normal(0,0.03))
		ys = [row.age, row.eeg_age]
		lines!(ax, xs, ys; color=:gray60, linestyle=:dash, linewidth=0.5)
		scatter!(ax, xs, ys; color=[(c, 0.3) for c in clrs], strokewidth=0.5, markersize=fsea_marker_size)
	end
end
	

#-

colsize!(figure3.layout, 2, Relative(4/5))
linkyaxes!(ax_future_violins...)


for (i, tp) in enumerate(ftps)
    @info "$tp"
    for (j, feat) in enumerate(filter(contains("latency"), eeg_features))
        gdf = groupby(futfsea_df, "eeg_feature")
        featstr = replace(feat, "peak_latency_"=>"", "_corrected"=>"c")
        @warn "$featstr"

		df =  alllms[(; eeg_feature=feat, timepoint=tp)]
		allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        for row in eachrow(subdf)
			yidx = df[!, row.geneset]
			xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
			ys = df.z[yidx]
			xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
            ymed = median(ys)
            
            cidx = row.q₀ > 0.2  ? 4 :       # not significant
                   row.q₀ < 0.01 ? 1 :       # quite significant
                   row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? colors_sig[cidx] : colors_sig[8 - cidx]
            violin!(ax_futlat[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_futlat[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
			lines!(ax_futlat[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end

    for (j, feat) in enumerate(filter(contains("amp"), eeg_features))
        gdf = groupby(futfsea_df, "eeg_feature")
        featstr = replace(feat, "peak_amp"=>"", "_corrected"=>"c")
        @warn "$featstr"

        df =  alllms[(; eeg_feature=feat, timepoint=tp)]
        allymed = median(df.z)

        subdf = subset(gdf[(; eeg_feature=feat)], "timepoint" => ByRow(==(tp)))
        # lines!(ax_futlat, [allymed, allymed], [0.5, size(subdf, 1)+0.5]; color=:gray, linestyle=:dash) 
        for row in eachrow(subdf)
            # haskey(gsidx, row.geneset) || continue
            yidx = df[!, row.geneset]
            xpos = gsidx[row.geneset] * gs_interval + (2 - i) * tp_interval
            ys = df.z[yidx]
            xs = rand(Normal(0.0, tp_interval / 8), length(ys)) .+ xpos
            ymed = median(ys)

            cidx = row.q₀ > 0.2  ? 4 :       # not significant
            row.q₀ < 0.01 ? 1 :       # quite significant
            row.q₀ < 0.1  ? 2 : 3 # somewhat significant / not very significant
            c = ymed < allymed ? colors_sig[cidx] : colors_sig[8 - cidx]
            violin!(ax_futamp[j][1], fill(xpos, length(ys)), ys; width=tp_interval*1.5, color=(c, 0.4), orientation=:horizontal)
            scatter!(ax_futamp[j][1], ys, xs; color = c, strokewidth=0.5, markersize = fsea_marker_size)    
            lines!(ax_futamp[j][1], [ymed, ymed], [xpos - tp_interval/2, xpos + tp_interval/2]; color = c)
        end
    end
    save("/home/kevin/Downloads/figure3_$tp.png", figure3)
end


#-

