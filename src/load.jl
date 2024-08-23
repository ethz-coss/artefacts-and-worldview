using DataFrames, DelimitedFiles
include("utils.jl")

# loading results methods
function processepsilon(path, e)
    res = readdlm(epsilonfilename(path, e), ',', Float64, '\n', header=true)
    df = DataFrame(res[1], vec(res[2]))
    df[!, :epsilon] .= e
    df.row_num = 1:nrow(df)
    stack(df, Not([:epsilon, :row_num]), view=true) 
end

function load_results(; epsilons, path)
    df = vcat(map(e->processepsilon(path, e), epsilons)...);
    df.type = last.(split.(df.variable));
    map!(i->chopsuffix(i, r" (alpha|beta|converged)$"), df.variable, df.variable);
    rename!(df, (:variable => :condition))
    df.type .= Symbol.(df.type)
    df[df.type .!= :converged, :], unstack(unique(df[df.type.==:converged, Not([:row_num, :type])]), :epsilon, :value)
    #df = unstack(df, :type, :value);
    #df[!, [:alpha, :beta]] = disallowmissing(df[!, [:alpha, :beta]])
end