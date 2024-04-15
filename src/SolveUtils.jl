
module SolveUtils
    using StatsBase, Distances, LinearAlgebra, SparseArrays
    export dist, F!, F, run

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
        x3::Array{T,1} = zeros(T, size(x,1))
    ) where {T<:Real, S<:Function}
        xprev .= x
        rows = rowvals(R)
        R_artefacts .= f.(pairwise(dist, artefacts, @view(x[:,1])))
        for (c, b) in enumerate(sum(R_artefacts, dims=1))
            rg = nzrange(R, c)
            rv = @view rows[rg]
            vs = map(i->f(dist(@view(xprev[rows[i], :]), @view(xprev[c, :]))), rg)
            s = sum(vs)
            x3[c] = cutoff((vs ./ s)' * @view(xprev[rv, 1]))
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
        x3::Array{T,1} = zeros(T, size(x,1))
    ) where {T<:Real, S<:Function}
        xprev .= x
        rows = rowvals(R)
        for c in 1:size(R,1)
            rg = nzrange(R, c)
            rv = @view rows[rg]
            vs = map(i->f(dist(@view(xprev[rows[i], :]), @view(xprev[c, :]))), rg)
            avg = (vs ./ sum(vs))'
            x3[c] = cutoff(avg * @view(xprev[rv, 1]))
            x[c,2] = cutoff(avg * (λ .* @view(xprev[rv, 2]) + (1-λ) .* @view(xprev[rv, 1])))
            x[c,1] = cutoff(avg * @view(xprev[rv, 1]))
        end
    end


    function Fconverged!(x, args...; xprev=zeros(size(x)), tol=1e-03, kwargs...)
        F!(x, args...; xprev, kwargs...)
        all(isapprox.(x, xprev; rtol=tol, atol=tol))
    end

    function run!(args...; iteration_steps = 50, check_every = 1, tol = 1e-05, kwargs...)
        i = 1
        while i <= iteration_steps
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
        x3i = zeros(N)
        xprev=zeros(size(xi))

        function process(n, a)
            @show n
            ia = 0
            xa = copy(xi)
            x3 = copy(x3i)
            if isempty(a)
                ia = run!(xa, R; λ, f, xprev, x3, kwargs...)
            else
                Ra = zeros((length(a), N))
                ia = run!(xa, R, Ra, a; λ, f, xprev, x3, kwargs...)
            end
            ia => (xa, x3)
        end
        Dict("initial" => (0 => (xi, x3i)), Tuple(n=>process(n, a) for (n, a) in conditions)...)
    end
    
end