function plotbugs!(ax, datadf, unirefs)

end

function topbugs(df, n=5; groupcol="genus")
    gdf = sort(combine(groupby(df, groupcol), "abundance" => sum => "abundance"), "abundance"; rev=true)
    notunclass = subset(gdf, groupcol => ByRow(!=("unclassified")))
    return size(notunclass, 1) > n ? notunclass[1:n, groupcol] : notunclass[!, groupcol]
end


function grouptop(df, n=5; groupcol="genus")
    keepbugs = topbugs(df, n; groupcol)

    combine(groupby(df, "sample"), AsTable([groupcol, "abundance"]) => (nt -> begin
        idx = findall(bug -> bug ∈ keepbugs, nt[Symbol(groupcol)])
        (; Symbol(groupcol) => [nt[Symbol(groupcol)][idx]; ["other"]],
            :abundance => [nt.abundance[idx]; [sum(nt.abundance[Not(idx)])]]
        )
    end) => [groupcol, "abundance"])
end
