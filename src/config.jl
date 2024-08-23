using Configurations, Random

@option mutable struct Config
    seed::Int
    rng::Maybe{<:Random.AbstractRNG}

    simpath::String
    figpath::String
    tabpath::String
    graphfilepath::String

    N::Int
    lambda::Float64
    conditionorder::Array{String}
    conditions::Dict{String, Vector{Float64}}
    epsilons::Array{Float64}

    t_max::Int = 12
    check_every::Int = 50
    tol::Float64 = 1e-06
    
end

Configurations.from_dict(::Type{Config}, ::Type{Array{Float64}}, s) = range(; [(Symbol(k)) => v for (k,v) in s]...)
config = from_toml(Config, "config.toml")
config.rng = Random.Xoshiro(config.seed)