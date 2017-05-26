import ArgParse, Gurobi, DataFrames
using ArgCheck

const arg_parse_settings = ArgParse.ArgParseSettings(prog="var", fromfile_prefix_chars='@',
                                                     description="Calculates the steady state as function of ξ.")

@ArgParse.add_arg_table arg_parse_settings begin
    "--sto"
        help = "stoichiometric matrix"
        arg_type = String
        required = true
    "--mets"
        help = "metabolites"
        arg_type = String
        required = true
    "--rxns"
        help = "reactions"
        arg_type = String
        required = true
    "--A"
        help = "biomass synthesis cost"
        arg_type = Float64
        default = 0.
        range_tester = A -> 0 ≤ A < Inf
    "--c"
        help = "media"
        arg_type = String
        required = true
    "--xi"
        help = "File with a list of ξ values to use. Execution stops at the first infeasible value of ξ found."
        arg_type = String
        required = true
    "--out"
        help = "output file"
        arg_type = String
        required = true
    "--prefix"
        help = "all paths arguments (--sto, --mets, --rxns, --c, --xi, --out) will be prefixed with this string."
        arg_type = String
        default = ""
end

const args = ArgParse.parse_args(arg_parse_settings)

set_multiobj_n!(gurobi_model::Gurobi.Model, n::Int) = Gurobi.set_intattr!(gurobi_model, "NumObj", n)

function set_multiobj!(gurobi_model::Gurobi.Model, i::Int, c, priority, weight)
    @assert length(c) == Gurobi.num_vars(gurobi_model)
    oldi = Gurobi.get_int_param(gurobi_model, "ObjNumber")
    Gurobi.set_int_param!(gurobi_model, "ObjNumber", i)
    Gurobi.set_dblattrarray!(gurobi_model, "ObjN", 1, Gurobi.num_vars(gurobi_model), c)
    Gurobi.set_intattr!(gurobi_model, "ObjNPriority", priority)
    Gurobi.set_dblattr!(gurobi_model, "ObjNWeight", weight)
    Gurobi.set_int_param!(gurobi_model, "ObjNumber", oldi)
end

function set_ub!(gurobi_model::Gurobi.Model, u::Vector{Float64})
    @assert length(u) == Gurobi.num_vars(gurobi_model)
    Gurobi.set_dblattrarray!(gurobi_model, "UB", 1, Gurobi.num_vars(gurobi_model), u)
end

function set_rhs!(gurobi_model::Gurobi.Model, rhs::Vector{Float64})
    @assert Gurobi.num_constrs(gurobi_model) == length(rhs)
    Gurobi.set_dblattrarray!(gurobi_model, "RHS", 1, Gurobi.num_constrs(gurobi_model), rhs)
end

function GurobiEnv()
    env = Gurobi.Env()
    # suppress Gurobi logging output
    Gurobi.setparam!(env, "OutputFlag", 0)
    Gurobi.setparam!(env, "Quad", 1)
    Gurobi.setparam!(env, "Method", 4)
    #Gurobi.setparam!(env, "MultiObjMethod", 1)
    Gurobi.setparam!(env, "FeasibilityTol", 1e-9)
    Gurobi.setparam!(env, "OptimalityTol", 1e-9)
    Gurobi.setparam!(env, "MarkowitzTol", 0.9)
    #Gurobi.setparam!(env, "LogToConsole", 0)
    #Gurobi.setparam!(env, "LogFile", path_prefix * "gurobi.log"))
    return env
end

function solve!(gurobi_model::Gurobi.Model)::Bool
    Gurobi.update_model!(gurobi_model)
    Gurobi.optimize(gurobi_model)
    return Gurobi.get_status(gurobi_model) == :optimal
end

function U(V, c, ξ)::Float64
    @argcheck V ≥ 0 && c ≥ 0 && 0 ≤ ξ < Inf
    return c == 0 ? 0. : min(V, c / ξ)
end


function main()
    path_prefix = strip(args["prefix"])

    ξlist = readdlm(path_prefix * strip(args["xi"]))
    @assert all(0 .≤ ξlist .< Inf)

    sto = DataFrames.readtable(path_prefix * strip(args["sto"]); allowcomments=true, separator=' ')
    S = sparse(sto[:i], sto[:k], sto[:s])
    m, n = Base.size(S)
    @assert m < n

    mets = DataFrames.readtable(path_prefix * strip(args["mets"]); allowcomments=true, separator=' ')
    rxns = DataFrames.readtable(path_prefix * strip(args["rxns"]); allowcomments=true, separator=' ')
    @assert Base.size(mets, 1) == m && Base.size(rxns, 1) == n

    @assert all(rxns[:lb] .≤ 0 .≤ rxns[:ub])
    @assert all(rxns[:an] .≥ 0) && all(rxns[:ap] .≥ 0) && all(mets[:bp] .≥ 0) && all(mets[:bn] .≥ 0)
    @assert all(mets[:L] .≤ 0 .≤ mets[:V])
    @assert all(-Inf .< mets[:e] .< Inf) && all(-Inf .< mets[:y] .< Inf)

    c = DataFrames.readtable(path_prefix * strip(args["c"]); allowcomments=true, separator=' ')
    @assert all(c[:mM] .≥ 0)
    @assert all(c[:id] == mets[:id])

    x = zeros(2n + 2m + 2)   # [r+, r-, u+, u-, μ, ϕ]
    upper::Vector{Float64} = [rxns[:ub]; -rxns[:lb]; mets[:V]; -mets[:L]; Inf; 1.]     # upper bounds
    @assert length(upper) == length(x) == 2m + 2n + 2

    # equality constrains
    W = sparse([ S            -S              speye(m)          -speye(m)     -mets[:y]       spzeros(m)
                 rxns[:ap]'    rxns[:an]'     mets[:bp]'        mets[:bn]'    args["A"]       -1 ])

    @assert Base.size(W) == (m + 1, length(x))

    # the problem has the form:
    # max μ s.t. W⋅x = [e,0] and 0 ≤ x ≤ upper.
    # with the secondary objective: min ϕv

    LP = Gurobi.gurobi_model(GurobiEnv(), sense = :maximize,
                             Aeq = W, beq = [collect(mets[:e]); 0],
                             lb = zeros(x), ub = upper, f = zeros(x))

    set_multiobj_n!(LP, 2)
    set_multiobj!(LP, 0, [zeros(2n + 2m); 1;  0], 10, 1.0)  # max μ
    set_multiobj!(LP, 1, [zeros(2n + 2m); 0; -1], 5,  1.0)  # min ϕ

    r = zeros(n)
    u = zeros(m)
    s = zeros(m)

    xout = open(path_prefix * strip(args["out"]) * ".x", "w")
    uout = open(path_prefix * strip(args["out"]) * ".u", "w")
    rout = open(path_prefix * strip(args["out"]) * ".r", "w")
    sout = open(path_prefix * strip(args["out"]) * ".s", "w")

    write(xout, "xi\tmu\tphi\n")
    write(uout, "xi"); for met in mets[:id]; write(uout, "\t" * met); end; write(uout, "\n")
    write(rout, "xi"); for rxn in rxns[:id]; write(rout, "\t" * rxn); end; write(rout, "\n")
    write(sout, "xi"); for met in mets[:id]; write(sout, "\t" * met); end; write(sout, "\n")

    for ξ in ξlist
        # update uptake bounds
        view(upper, 2n + (1:m)) .= U.(mets[:V], c[:mM], ξ)
        set_ub!(LP, upper)

        solve!(LP) || return
        x .= Gurobi.get_solution(LP)
        r .= x[1:n] .- x[n + (1:n)]
        u .= x[2n + (1:m)] .- x[2n + m + (1:m)]
        s .= c[:mM] .- u * ξ
        μ = x[2n + 2m + 1]
        ϕ = x[2n + 2m + 2]

        write(xout, "$ξ\t$μ\t$ϕ\n")
        @printf "ξ=%-7.13f μ=%-7.13f ϕ=%-7.13f\n" ξ μ ϕ

        write(uout, "$ξ"); for i = 1:m; write(uout, "\t$(u[i])"); end; write(uout, "\n")
        write(rout, "$ξ"); for k = 1:n; write(rout, "\t$(r[k])"); end; write(rout, "\n")
        write(sout, "$ξ"); for i = 1:m; write(sout, "\t$(s[i])"); end; write(sout, "\n")
    end

    close(xout)
    close(uout)
    close(rout)
    close(sout)
end


main()
