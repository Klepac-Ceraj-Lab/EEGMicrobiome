# https://github.com/JuliaStats/MultivariateStats.jl/pull/162
function loadings(M::MultivariateStats.MDS)
    ev = eigvals(M)
    return ev' .* projection(M)[:, 1:length(ev)]
end

# https://github.com/JuliaStats/MultivariateStats.jl/pull/162
function loadings(M::MultivariateStats.MDS, dim)
    l = loadings(M)
    return l[:, dim]
end

varexplained(M::MultivariateStats.MDS) = eigvals(M) ./ sum(eigvals(M))

mdsaxis(M::MultivariateStats.MDS, dim::Int) = "MDS$dim ($(round(varexplained(M)[dim] * 100, digits=2))%)"

function plot_pcoa!(ax, M::MultivariateStats.MDS; dims=(1,2), kwargs...)
    ax.xlabel = get(kwargs, :xlabel, mdsaxis(M, dims[1]))
    ax.ylabel = get(kwargs, :ylabel, mdsaxis(M, dims[2]))
    
    scatter!(ax, loadings(M, dims[1]), loadings(M, dims[2]); kwargs...)
end

