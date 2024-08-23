using DataFrames, DelimitedFiles, SparseArrays, Printf, StatsBase, Distances, LinearAlgebra
include("config.jl")
include("utils.jl")

dist = euclidean
cutoff(x, lower=0, upper=1) = max(min(x, upper), lower)

function F!(
    x::Array{T, 2}, 
    R::SparseMatrixCSC{<:Any, Int},
    R_artefacts::Array{T, 2}, 
    artefacts::AbstractArray{T, 1};
    xprev::Array{T, 2}, 
    f::S,
    λ::Real,
) where {T<:Real, S<:Function}
    xprev .= x
    rows = rowvals(R)
    R_artefacts .= f.(pairwise(dist, artefacts, @view(x[:,1])))
    for (c, b) in enumerate(sum(R_artefacts, dims=1))
        rg = nzrange(R, c)
        rv = @view rows[rg]
        vs = map(i->f(dist(@view(xprev[rows[i], :]), @view(xprev[c, :]))), rg)
        s = sum(vs)
        x[c,2] = cutoff((vs ./ s)' * (λ .* @view(xprev[rv, 2]) + (1-λ) .* @view(xprev[rv, 1])))
        s += b
        x[c,1] = cutoff((vs ./ s)' * @view(xprev[rv, 1]) + (@view(R_artefacts[:, c]) ./s)' * artefacts)
    end
end

function F!(
    x::Array{T, 2}, 
    R::SparseMatrixCSC{<:Any, Int};
    xprev::Array{T, 2}, 
    f::S,
    λ::Real,
) where {T<:Real, S<:Function}
    xprev .= x
    rows = rowvals(R)
    for c in 1:size(R,1)
        rg = nzrange(R, c)
        rv = @view rows[rg]
        vs = map(i->f(dist(@view(xprev[rows[i], :]), @view(xprev[c, :]))), rg)
        avg = (vs ./ sum(vs))'
        x[c,2] = cutoff(avg * (λ .* @view(xprev[rv, 2]) + (1-λ) .* @view(xprev[rv, 1])))
        x[c,1] = cutoff(avg * @view(xprev[rv, 1]))
    end
end


function Fconverged!(x, args...; xprev=zeros(size(x)), tol=1e-03, kwargs...)
    F!(x, args...; xprev, kwargs...)
    all(isapprox.(x, xprev; rtol=tol, atol=tol))
end


function run!(args...; t_max = 50, check_every = 1, tol = 1e-05, kwargs...)
    i = 1
    while i <= t_max
        if i % check_every == 0 && Fconverged!(args...; tol, kwargs...)
            return i
        end
        F!(args...; kwargs...)
        i+= 1
    end
    Fconverged!(args...; tol, kwargs...) ? i : -1
end


function run(rng; R, N, λ, f, conditions, kwargs...)
    
    xi = rand(rng, N, 2)
    xprev=zeros(size(xi))

    function process(n, a)
        ia = 0
        xa = copy(xi)
        if isempty(a)
            ia = run!(xa, R; λ, f, xprev, kwargs...)
        else
            Ra = zeros((length(a), N))
            ia = run!(xa, R, Ra, a; λ, f, xprev, kwargs...)
        end
        ia => xa
    end
    Dict("initial" => (0 => xi), Tuple(n=>process(n, a) for (n, a) in conditions)...)
end

function social_nw_R(graphfilename, N; T=Bool, directed=true)
    df = readdlm(graphfilename, ',', Float64, '\n')
    if directed
        sparse(vcat(df[:, 1],1:N), vcat(df[:, 2],1:N), zeros(T, size(df, 1)+N), N, N)
    else 
        sparse(vcat(df[:, 1],df[:, 2],1:N), vcat(df[:, 2],df[:, 1],1:N), zeros(T, 2*size(df, 1)+N), N, N)
    end
end

function write_artefact_epsilons(R, N; rng,
    λ,
    path,
    epsfunc=(i::Real, e::Real) -> Int(i ≤ e), 
    conditions,
    epsilons,
    kwargs...
)
    
    colnames = ["$a $c" for a in ("initial",  keys(conditions)...) for c in ["alpha", "beta", "converged"]]
    df = DataFrame([i => Vector{Float64}(undef, N) for i in colnames], copycols=false)
    for e in epsilons
        x_artefacts = run(rng; R, N, λ, f= i-> epsfunc(i,e), conditions, kwargs...)
        for (c, (i, x)) in x_artefacts
            @view(df[:, ["$(c) alpha", "$(c) beta"]]) .= x
            @view(df[:, ["$(c) converged"]]) .= i
        end
        writedlm(epsilonfilename(path, e), Iterators.flatten(([names(df)], eachrow(df))), ',')
    end
end

if !ispath(config.simpath); mkdir(config.simpath); end
R = social_nw_R(config.graphfilepath, config.N; T=Bool, directed=true)
write_artefact_epsilons(
    R, config.N;
    λ = config.lambda, 
    path=config.simpath,
    rng=config.rng, 
    epsilons=config.epsilons, 
    conditions=config.conditions, 
    t_max=config.t_max, 
    check_every=config.check_every, 
    tol=config.tol
)