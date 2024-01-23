
function runlms(indf, outfile, respcol, featurecols;
    age_col = "age",
    additional_cols = Symbol[],
    formula = term(:func) ~ term(respcol) + term(age_col) + term(:age_diff) + term(:n_segments),
    prevalence_filter = 0.05
    )
    @debug "Respcol: $respcol"
    lmresults = DataFrame(ThreadsX.map(featurecols) do feature
    # @debug "Feature: $feature"

    # ab = collect(indf[!, feature] .+ (minimum(indf[over0, feature])) / 2) # add half-minimum non-zerovalue
    df = select(indf, respcol, age_col, "age_diff", "n_segments", additional_cols...)
    over0 = indf[!, feature] .> 0
    default_ret = (; feature, Name = respcol, coef = NaN, std_err = NaN, z = NaN, pvalue = NaN, lower_95 = NaN, upper_95 = NaN, qvalue = NaN)

    (sum(skipmissing(over0)) / length(over0)) > prevalence_filter || return default_ret

    df.func = over0
    # @debug "DataFrame: $df"

    try
    mod = glm(formula, df, Binomial(), LogitLink(); dropcollinear=false)
    ct = DataFrame(coeftable(mod))
    ct.feature .= feature
    rename!(ct, "Pr(>|z|)"=>"pvalue", "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
    select!(ct, Cols(:feature, :Name, :))
    return NamedTuple(only(filter(row-> row.Name == respcol, eachrow(ct))))    
    catch e
    @warn "hit $e for $feature"
    return default_ret
    end
    end)

    subset!(lmresults, "pvalue"=> ByRow(!isnan))
    DataFrames.transform!(lmresults, :pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(lmresults, :qvalue)

    CSV.write(outfile, lmresults)
    lmresults
end
