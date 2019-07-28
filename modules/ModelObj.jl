module ModelObj

using Parameters
using CSV
using JLD
using DataFrames

bpr_inputs_dict = Dict{Symbol, Real}(:vtmax => 1.2,
                                     :vtN => 15,
                                     :ttmN => 10, 
                                     :vbhlmin => .6,
                                     :vbhlmax => 1.4,
                                     :vbhlN => 11,
                                     :vmax => 1.5,
                                     :vN => 10^3,
                                     :uN => 10^3,
                                     :vtN_ref => 600,
                                     :ttmN_ref => 450,
                                     :vbhlN_ref => 550)

obj_params = [:mu_b, :m, :c, :p, :vbl, :vbh]
firm_params = [:V0, :alpha, :pi, :r, :gross_delta, :iota,
               :xi, :kappa, :lambda, :sigmal, :sigmah]


# * Model Structs ##############################################
# include("model_structs.jl")
@with_kw mutable struct KStruct
    mu_b::Float64
    m::Float64
    c::Float64
    p::Float64
    vbl::Float64
end


@with_kw mutable struct FirmParams
    V0::Float64
    alpha::Float64
    pi::Float64
    
    r::Float64
    gross_delta::Float64
    iota::Float64
    
    xi::Float64
    kappa::Float64
    
    lambda::Float64
    sigmal::Float64
    sigmah::Float64
end


# Bond Pricing Inputs
@with_kw mutable struct BPrInputs
    # Grid Bounds & Lengths
    vtmax::Float64
    vtN::Int64
    ttm_max::Float64
    ttmN::Int64
    vbhlmin::Float64
    vbhlmax::Float64
    vbhlN::Int64

    vmax::Float64
    vN::Int64
    uN::Int64

    vtN_ref::Int64
    ttmN_ref::Int64
    vbhlN_ref::Int64
end


# Bond Pricing Inputs
@with_kw mutable struct BPrFixedTTMInputs
    # Grid Bounds & Lengths
    vtmax::Float64
    vtN::Int64
    ttm::Float64
    vbhlmin::Float64
    vbhlmax::Float64
    vbhlN::Int64
end


# Bond Pricing Grids & Surfaces
@with_kw mutable struct BPrSurfs
    # Grids
    vtgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    ttmgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    vbhlgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    
    # Surfs
    f11_surf::Array{Float64,2}
    f12_surf::Array{Float64,3}
    f13_surf::Array{Float64,3}
    f21_surf::Array{Float64,2}
    f22_surf::Array{Float64,2}
end


# Bond Pricing Interpolated Functions
@with_kw mutable struct BPrInterpFuns
    f11
    f12
    f13
    f21
    f22
end


@with_kw mutable struct Firm
    mu_b::Float64
    m::Float64
    c::Float64
    p::Float64
    vbl::Float64
    vbh::Float64
    pm::FirmParams
    model::String
    bi::BPrInputs
    bs::BPrSurfs
    bf::BPrInterpFuns
    bit::BPrFixedTTMInputs
    bft::BPrInterpFuns
    optKS::KStruct
end


# * Set Functions ##############################################
# include("model_set_funs.jl")
function grid_creator(z0, z1, n)
    dz = (z1 - z0) / float(n)
    if z0^2 < 1e-6
        return dz, range(z0+dz, stop=z1, length=n)
    else
        return dz, range(z0, stop=z1, length=n)
    end
    # return dz, range(z0, stop=z1, length=n) 
end


function set_bpr_inputs_struct(bprf)
    return BPrInputs(bprf[:vtmax], 
                     bprf[:vtN],
                     bprf[:ttm_max],
                     bprf[:ttmN],
                     bprf[:vbhlmin],
                     bprf[:vbhlmax],
                     bprf[:vbhlN],
                     bprf[:vmax],
                     bprf[:vN], 
                     bprf[:uN], 
                     bprf[:vtN_ref],
                     bprf[:ttmN_ref],
                     bprf[:vbhlN_ref])
end


function set_bpr_surfs_struct(bprf)
    return BPrSurfs(bprf[:vtgrid],
                    bprf[:ttmgrid],
                    bprf[:vbhlgrid],
                    bprf[:f11_surf],
                    bprf[:f12_surf], 
                    bprf[:f13_surf], 
                    bprf[:f21_surf], 
                    bprf[:f22_surf])
end


function set_bpr_grids(svm)
    _, svm.bs.vtgrid = grid_creator(0.0, svm.bi.vtmax, svm.bi.vtN)
    _, svm.bs.ttmgrid = grid_creator(0.0, svm.bi.ttm_max, svm.bi.ttmN)
    _, svm.bs.vbhlgrid = grid_creator(svm.bi.vbhlmin, svm.bi.vbhlmax, svm.bi.vbhlN)

    return svm
end


function set_bpr_surfs(svm, bprf)
    svm.bs.f11_surf = bprf[:f11_surf]
    svm.bs.f12_surf = bprf[:f12_surf]
    svm.bs.f13_surf = bprf[:f13_surf]
    svm.bs.f21_surf = bprf[:f21_surf]
    svm.bs.f22_surf = bprf[:f22_surf]

    return svm
end


function set_bpr_funs_struct(svm, bprf)
    svm.bf.f11 = bprf[:f11]
    svm.bf.f12 = bprf[:f12]
    svm.bf.f13 = bprf[:f13]
    svm.bf.f21 = bprf[:f21]
    svm.bf.f22 = bprf[:f22]

    return svm
end


function set_opt_k_struct(mobj, df::DataFrame; m::Float64=NaN, tol::Float64=1e-6)

    if mobj.model == "svm"
        vars = fieldnames(ModelObj.FirmParams)
    elseif mobj.model == "cvm"
        vars = [x for x in fieldnames(ModelObj.FirmParams) if
                !(x in [:lambda, :sigmah])]
    else
        println("Please enter either a SVM or a CVM object. Exiting...")
        return
    end

    # Find Index
    LL = []
    for var in vars 
        append!(LL, [abs.(df[:, var] .- getfield(mobj.pm, var)) .< tol])
    end

    if !isnan(m)
        append!(LL, [abs.(df[:, :m] .- m) .< tol])
    end
    index = .&(LL...)


    # Set Capital Structure
    for var in [:mu_b, :m, :c, :p]
        setfield!(mobj.optKS, var, df[index, var][1])
    end
    mobj.optKS.vbl = df[index, :vb][1]
    
    return mobj
end


# * Load Functions #############################################
# include("model_load_funs.jl")
function load_firm_params_bpr_inputs(svm_input_path::String; file_name::String="_svm_inputs")
    df = CSV.read(string(svm_input_path, "/", file_name, ".csv"))

    for var in obj_params
        if !(var in names(df))
            df[!, var] = NaN
        end
    end
    
    # ####### Generate Firm Params Dictionary #######
    params_list = vcat(obj_params, firm_params)    
    params = Dict{Symbol, Any}(cn => df[1, cn] for cn in params_list if
                                !occursin(:Column, cn))
    
    # ####### Generate BPr Inputs Dictionary #######
    cols = [:vtN, :ttmN, :vbhlN, :vN, :uN]
    bpr_inputd = Dict{Symbol, Any}(cn => df[1, cn] for cn in names(df) if
                                   !occursin(:Column, cn) &
                                   !(cn in cols) & !(cn in params_list))

    # Grid Length Numbers Should be Interger
    for cn in cols
        bpr_inputd[cn] = trunc(Integer, df[1, cn])
    end

    return  params, bpr_inputd
end


function load_bpr_surfs(svm_input_path::String)
    fsurfs = [string("f", f, "_surf") for f in ["11", "12", "13", "21", "22"]]
    bprf = Dict()
    for f in fsurfs
        bprf[Symbol(f)] = load(string(svm_input_path, "/", f, ".jld"))[f]
    end

    return bprf
end


# * Extract Parameters and Get Object ##########################
function extract_param(svm, pname::Symbol)
    val = NaN
    
    svm_pos = findall(x -> x==pname, [x for x in fieldnames(typeof(svm))])
    svm_pm_pos = findall(x -> x==pname, [x for x in fieldnames(typeof(svm.pm))])
    
    if length(svm_pos) > 0
        # return getfield(svm, fieldnames(svm)[svm_pos[1]])
        return getfield(svm, fieldnames(typeof(svm))[svm_pos[1]])
    else
        # return getfield(svm.pm, fieldnames(svm.pm)[svm_pm_pos[1]])
        return getfield(svm.pm, fieldnames(typeof(svm.pm))[svm_pm_pos[1]])
    end
end

function get_obj_model(svm)
    return any([svm.model == "cvm",
                isnan(svm.pm.lambda),
                isnan(svm.pm.sigmah)]) ? "cvm" : "svm"
end

# * Model Constructors ########################################
# include("model_constructors.jl")
function firm_constructor(dict::Dict{Symbol,Float64};
                          bpr_inputs::Dict{Symbol, Real}=bpr_inputs_dict,
                          model::String="svm")
    for par in [:vbl, :vbh, :c, :p]
        if !(par in keys(dict))
            println(string("Setting initial ", par, " value to ", NaN))
            dict[par] = NaN
        end
    end

    if !(:mu_b in keys(dict)) | isnan(dict[:mu_b])
        mu_b = 1.0
        println(string("Setting measure of bonds to ", mu_b))
        dict[:mu_b] = mu_b
    end

    if model == "cvm"
        println("Constant Volatility Model: setting vbh to vbl, lambda to NaN")
        dict[:lambda] = NaN
        dict[:vbh] = dict[:vbl]
    end

    # Set Firm Parameters (dynamics, default, liquidity)
    firm_params =  FirmParams(dict[:V0], dict[:alpha], dict[:pi],
	                      dict[:r], dict[:gross_delta], dict[:iota],
	                      dict[:xi], dict[:kappa],
	                      dict[:lambda], dict[:sigmal], dict[:sigmah])

    # Bond Pricing Inputs
    ttm_max = dict[:m]
    bpr_surfs_inputs = BPrInputs(bpr_inputs[:vtmax], 
                                 bpr_inputs[:vtN],
                                 ttm_max,
                                 bpr_inputs[:ttmN],
                                 bpr_inputs[:vbhlmin],
                                 bpr_inputs[:vbhlmax],
                                 bpr_inputs[:vbhlN],
                                 bpr_inputs[:vmax],
                                 bpr_inputs[:vN], 
                                 bpr_inputs[:uN], 
                                 bpr_inputs[:vtN_ref],
                                 bpr_inputs[:ttmN_ref],
                                 bpr_inputs[:vbhlN_ref])


    # Bond Pricing Grids & Surfaces
    _, vtgrid = grid_creator(0.0, bpr_surfs_inputs.vtmax, bpr_surfs_inputs.vtN)
    _, ttmgrid = grid_creator(0.0, bpr_surfs_inputs.ttm_max, bpr_surfs_inputs.ttmN)
    _, vbhlgrid = grid_creator(bpr_surfs_inputs.vbhlmin, bpr_surfs_inputs.vbhlmax, bpr_surfs_inputs.vbhlN)
        
    mat1 = fill(NaN, 1, 1)
    mat2 = fill(NaN, 1, 1, 1)
        
    bpr_surfs = BPrSurfs(vtgrid, ttmgrid, vbhlgrid,
                     mat1, mat2, mat2, mat1, mat1)

    # Bond Pricing Interpolated Functions 
    bpr_interp_funs = BPrInterpFuns(nothing, nothing, nothing,
                                    nothing, nothing)

    bpr_fixed_ttm_inputs = BPrFixedTTMInputs(NaN, 0, NaN, NaN, NaN, 0)

    bpr_fixed_ttm_interp_funs = BPrInterpFuns(nothing, nothing, nothing,
                                nothing, nothing)

    # Optimal Capital Structs
    optKS = KStruct(NaN, NaN, NaN, NaN, NaN) 
  

    # Construct Firm Object
    return Firm(dict[:mu_b], dict[:m],
                dict[:c], dict[:p],
	    	dict[:vbl], dict[:vbh],
	        firm_params, model,
                bpr_surfs_inputs,
                bpr_surfs, 
                bpr_interp_funs,
                bpr_fixed_ttm_inputs,
                bpr_fixed_ttm_interp_funs,
                optKS)

end


function create_misrep_type(svmi, svmj)
    # Extract Type-i parameters
    ijdict = Dict{String, Float64}()
    for par in vcat(obj_params, firm_params)
        ijdict[string(par)] = extract_param(svmi, string(par))
    end

    # Set iota to zero
    ijdict[:iota] = .0

    # Copy capital structure from Type-j
    for par in [:mu_b, :c, :p]
        ijdict[par] = extract_param(svmj, par)
    end

    return firm_constructor(ijdict)
end


# * END MODULE
end
