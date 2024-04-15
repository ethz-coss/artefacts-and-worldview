module LoadUtils
    include("SolveUtils.jl")
    using DataFrames, CSV, SparseArrays

    export social_nw_R, load_adf, write_artefact_epsilons

    function social_nw_R(graphfilename; N=10000, T=Bool, directed=true)
        df = CSV.read(graphfilename, DataFrame, header=false)
        if directed
            sparse(vcat(df[:, 1],1:N), vcat(df[:, 2],1:N), zeros(T, size(df, 1)+N), N, N)
        else 
            sparse(vcat(df[:, 1],df[:, 2],1:N), vcat(df[:, 2],df[:, 1],1:N), zeros(T, 2*size(df, 1)+N), N, N)
        end
    end

    # writing simulation results
    function write_artefact_epsilons(R, N, λ; rng, fname,
        epsfunc=(i::Real, e::Real) -> Int(i ≤ e), 
        conditionorder,
        conditions,
        epsilons=.05:.01:.4,
        iteration_steps=50,
        check_every=2,
        tol=1e-05
    )
        colnames = ["$a $c" for a in ("initial",  conditionorder...) for c in ["alpha", "beta", "x3"]]
        df = DataFrame([i => Vector{Float64}(undef, N) for i in colnames], copycols=false)
        for e in epsilons
            x_artefacts = SolveUtils.run(rng; R, N, λ, f= i-> epsfunc(i,e), conditions, iteration_steps, check_every, tol)
            for (c, (_, (x, x3))) in x_artefacts
                @view(df[:, ["$(c) alpha", "$(c) beta"]]) .= x
                @view(df[:, "$(c) x3"]) .= x3
            end
            print("epsilon=$(e).csv converged after $(keys(x_artefacts) .=> first.(values(x_artefacts)))")
            CSV.write(joinpath(fname, "epsilon=$(e).csv"), df)
        end
    end

    # loading results methods
    function processepsilon(ei, path)
        df = CSV.read(joinpath(path, "epsilon=$(ei).csv"), DataFrame) 
        df[!, :epsilon] .= ei
        df.row_num = 1:nrow(df)
        stack(df, Not(:epsilon, :row_num), view=true) 
    end

    function load_adf(args...; e =.05:.01:.4, path, kwargs...)
        adf = vcat(map(ei->processepsilon(ei, path), e)...);
        adf.type = last.(split.(adf.variable));
        map!(i->chopsuffix(i, r" (alpha|beta|x3)$"), adf.variable, adf.variable);
        adf = unstack(adf, :type, :value);
        adf[!, [:alpha, :beta]] = disallowmissing(adf[!, [:alpha, :beta]])
        rename!(adf, (:variable => :condition))
        adf
    end
end