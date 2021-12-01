module_path = "/home/artur/BondPricing/bond-model/modules/"
push!(LOAD_PATH, module_path)
# modnames = ["ModelObj", "AnalyticFunctions", "BondPrInterp", "EqFinDiff"]
# for modl in modnames
#    if !(joinpath(module_path, modl) in LOAD_PATH)
#        push!(LOAD_PATH, joinpath(module_path, modl))
#    end
# end


module Batch

using Parameters
using Printf
using Dierckx
using Distributed
using DataFrames
using JLD
using CSV
using Dates
using LinearAlgebra
using DSP
using Distributed


using ModelObj: firm_constructor
using AnalyticFunctions: get_cvm_vb, get_agg_c,
                         get_agg_p, get_rdisc,
                         get_param,
                         get_leverage, get_mbr
using BondPrInterp: get_cvm_debt_price, get_svm_debt_price,
                    bpr_surfs, bpr_interp, 
                    get_bond_yield, get_bond_spread
using EqFinDiff: get_cvm_eq, eq_fd

# * Batch Inputs
# include("_batch_inputs.jl")
# Batch Inputs
# Directories:
main_dir = "bond-model"
res_dir = "Results"
mat_dir_prefix = "m_"
coupon_dir_prefix = "c_"

batch_file_name = "batch_obj"

common_params = Dict{Symbol, Float64}(:V0 => 100.,
                                      :alpha => .6,
                                      :pi => .27,
                                      :r => .08)


_params_order = [:mu_b, :m, :xi, :kappa, :gross_delta,
                 :lambda, :iota, :sigmal, :sigmah]

# #######################################################
# ################# Combination Folders #################
# #######################################################
# Notice that mu_b is set to 1.
comb_folder_dict = Dict{Symbol,Array{String,1}}(:mu_b => ["mub_", "%.2f"],
                                                :m => ["m_", "%.2f"],
                                                :c => ["c_", "%.2f"],
                                                :xi => ["__xi_", "%.2f"],
                                                :kappa => ["__kappa_", "%.2f", "_bp"],
                                                :gross_delta => ["__gross_delta_", "%.2f", "_bp"],
                                                :iota => ["__iota_", "%.2f", "_bp"],
                                                :lambda => ["__lambda_", "%.3f"],
                                                :sigmal => ["__sigmal_", "%.3f"],
                                                :sigmah => ["__sigmah_", "%.3f"])


# #######################################################


# Define Parameter Combinations:
# #######################################################
# ######################### CVM #########################
# #######################################################
# Filename
cvm_sol_fn_prefix = "cvm_res"

# Parameters
cvm_param_values_dict = Dict{Symbol,Array{Float64,1}}(:mu_b => [1.],
                                                      :m => [1.], #, 3., 5., 7., 10.],
                                                      :xi => [1.],
                                                      :kappa => [k * 1e-4 for k in [25, 50, 100, 150]],
                                                      :gross_delta => [.02],
                                                      :lambda => [NaN],
                                                      :iota => [i * 1e-4 for i in [0.0, .5, 1., 1.5, 2.0, 2.5,
                                                                                   5, 7.5, 10, 15, 20, 25, 50]],
                                                      :sigmal => [.15, .2, .25, .30],
                                                      :sigmah => [NaN])

# Coupon Values Grid:
cvm_c_step = .5
cvm_coupon_grid = vcat([.25, .5], range(1, stop=25 + cvm_c_step, step=cvm_c_step))

# #######################################################
# ######################### SVM #########################
# #######################################################
# Julia Paths & Files


# FileNames
# svm_partI_fn_prefix = "svm_data_c"
# svm_partII_fn_prefix = "svm_res_c"
debt_at_par_cp_fn_prefix = string("debt_at_par_", coupon_dir_prefix)
eq_fd_cp_fn_prefix = string("eq_fd_", coupon_dir_prefix)




# Parameters
# svm_param_values_dict = Dict{Symbol,Array{Float64,1}}(:mu_b => [1.],
#                                                       :m => [1.],     #, 3., 5., 7., 10.],
#                                                       :xi => [1.],
#                                                       :kappa => [k * 1e-4 for k in [10, 25, 30, 40, 50]],
#                                                       :gross_delta => [.02],
#                                                       :lambda => [.1, .2, .3, .5, .75],
#                                                       :iota => [.0],
#                                                       :sigmal => [.15],
#                                                       :sigmah => [.175, .2, .225, .25, .275, .30, .325, .35])

svm_param_values_dict = Dict{Symbol,Array{Float64,1}}(:mu_b => [1.],
                                                    :m => [1., 3., 5., 7., 10.], # [.1]
                                                    :xi => [1.],
                                                    :kappa => [k * 1e-4 for k in [10, 25, 30, 40, 50]],
                                                    :gross_delta => [.02],
                                                    :lambda => [0.05, .1, .15, .2, .25, .3, .5, .75],
                                                    :iota => [.0],
                                                    :sigmal => [.15],
                                                    :sigmah => [.16, .175, .19, .2, .21, .225, .275, .30, .325, .35])



# Coupon Values Grid:
svm_c_step = .5
svm_coupon_grid = Array(range(svm_c_step, stop=11., step=svm_c_step))


# ### Set Tolerance Levels for Filtering p values###
# p_interp_fun -> eq_deriv_tol/debt_tol  -> p_tol_search
# eq_deriv_tol = vcat([1, 5] .* 1e-6,
#                     kron([1e-5, 1e-4, 1e-3, 1e-2],[1, 2.5, 5, 7.5]),
#                     [1, 2.5, 5] .* 1e-2)
# debt_tol = kron([1, 5], [10^i for i in range(0, stop=3)] .* 1e-4)
eq_deriv_tol1 = kron([1, 5], [1e-6, 1e-5, 1e-4])
debt_tol1 =  kron([1, 5], [10^i for i in [0, 1.]] .* 1e-4)
eq_deriv_tol2 = vcat(kron([1e-3, 1e-2],[1, 2.5, 5, 7.5]),
                    [1, 2.5, 5] .* 1e-2)
debt_tol2 = kron([1, 2.5, 5.], [10^i for i in range(0, stop=2)] .* 1e-4)

function toldf_fun(debt_tol::Array{Float64,1}, eq_deriv_tol::Array{Float64,1})
    # Form all possible combinations and sort them so that
    # debt is relaxed first.
    value_lists = [sort(debt_tol), sort(eq_deriv_tol)]
    cols = [:debt_diff, :eq_deriv]

    A = vcat([(i, j) for i=debt_tol, j=eq_deriv_tol]...)
    return sort(DataFrame(map(idx -> getindex.(A, idx),
                              eachindex(first(A))), cols),
                cols)
end
toldf = sort(vcat(toldf_fun(debt_tol1, eq_deriv_tol1),
                  toldf_fun(debt_tol2, eq_deriv_tol2)),
             [:eq_deriv, :debt_diff])



# Filter (P, VB) values -> discard candidates for which
# |(Debt - P)/P| > tol:
pvb_tol_vec = Array([1, 10, 10^2, 10^3]) * 10^-6



# DataFrame Columns
main_params = [:gross_delta, :delta, :iota, :kappa, :lambda, :sigmah, :sigmal]
k_struct_params = [:mu_b, :m, :c, :p, :vb]
fixed_params = [:V0, :r, :alpha , :pi, :xi]
debt_vars = [:debt_diff,
             :abs_debt_diff,
             :debt_per_diff,
             :abs_debt_per_diff]
equity_vars = [:eq_deriv,
               :eq_deriv_min_val,
               :eq_min_val,
               :eq_negative,
               :eq_vb]
# share_values = [:debt, :equity, :firm_value, :leverage, :MBR]
share_values = [:debt, :equity, :firm_value, :leverage, :MBR, :yield, :yield_spd]
dfcols = vcat(main_params, k_struct_params, debt_vars[1],
              equity_vars, debt_vars[2:end], share_values, fixed_params)



# #######################################################
# ####################### RESULTS #######################
# #######################################################
soldf_name = "soldf"
optdf_name = "optdf"
opt_k_struct_df_name="opt_k_structs"
#comb_opt_k_struct_df_coltypes=vcat(String, fill(Float64, 33), Bool, Float64)
opt_k_struct_df_coltypes=vcat(Int64, Float64, Int64, String, fill(Float64, 32), Bool, fill(Float64, 3))
cvm_opt_k_struct_df_coltypes=vcat([Int64, Float64, Int64, String],
                                  fill(Float64, 32))
# #######################################################


# * Batch Structs
# include("_batch_structs.jl")
@with_kw mutable struct DirFileNames
    # Directories
    main_dir::String
    res_dir::String
    mat_dir_prefix::String
    coupon_dir_prefix::String

    # File Names
    debt_at_par_cp_fn_prefix::String
    eq_fd_cp_fn_prefix::String
end


@with_kw mutable struct ModelInputs
    _svm_dict::Dict{Symbol, Float64}

    # Directories
    batch_res_dir::String
    maturity_dir::String
    comb_res_dir::String

    # Paths
    main_dir_path::String
    batch_res_path::String
    maturity_path::String
    comb_res_path::String
end


@with_kw mutable struct BatchDFColumns
    main_params::Array{Symbol,1}
    k_struct_params::Array{Symbol,1}
    fixed_params::Array{Symbol,1}
    debt_vars::Array{Symbol,1}
    equity_vars::Array{Symbol,1}
    share_values::Array{Symbol,1}
    dfcols::Array{Symbol,1}
end


@with_kw mutable struct BatchParams
    _param_values_dict::Dict{Symbol,Array{T,1} where T}
    _common_params::Dict{Symbol,Float64}
    _params_order::Array{Symbol,1}
    _params_combs::Array{Array{Float64,1},9}
    df::DataFrame
end



@with_kw mutable struct BatchStruct
    # Model Type
    model::String

    # Paths & File Names
    dfn::DirFileNames

    # DataFrame Columns
    dfc::BatchDFColumns

    # Model Object Inputs
    mi::ModelInputs

    # Parameters
    bp::BatchParams

    # Coupon
    # coupon_grid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
    #                           Base.TwicePrecision{Float64}}
    coupon_grid::Array{Float64,1}
    pvb_tol_vec::Array{Float64,1}
end


# * Batch Object
function BatchObj(;model::String="svm",
                  main_dir::String=main_dir,
                  res_dir::String=res_dir,
                  mat_dir_prefix::String=mat_dir_prefix,
                  coupon_dir_prefix::String=coupon_dir_prefix,
                  debt_at_par_cp_fn_prefix::String=debt_at_par_cp_fn_prefix,
                  eq_fd_cp_fn_prefix::String=eq_fd_cp_fn_prefix,
                  main_params::Array{Symbol,1}=main_params,
                  k_struct_params::Array{Symbol,1}=k_struct_params,
                  fixed_params::Array{Symbol,1}=fixed_params,
                  debt_vars::Array{Symbol,1}=debt_vars,
                  equity_vars::Array{Symbol,1}=equity_vars,
                  share_values::Array{Symbol,1}=share_values,
                  dfcols::Array{Symbol,1}=dfcols,
                  param_values_dict=svm_param_values_dict,
                  params_order::Array{Symbol,1}=_params_order,
                  coupon_grid::Array{Float64,1}=svm_coupon_grid,
                  pvb_tol_vec::Array{Float64,1}=pvb_tol_vec)

    if !(model in ["cvm", "svm"])
        println("Please choose model (cvm or svm.) Exiting...")
        return
    end

    if model == "cvm"
        debt_at_par_cp_fn_prefix = ""
        eq_fd_cp_fn_prefix = ""
        equity_vars = [:eq_min_val, :eq_vb]
        dfcols = vcat(main_params, k_struct_params, debt_vars[1],
                      equity_vars, debt_vars[2:end], share_values, fixed_params)
        param_values_dict = cvm_param_values_dict
        coupon_grid = cvm_coupon_grid
        pvb_tol_vec = Array{Float64,1}()
    end

    dfn = DirFileNames(main_dir,
                       res_dir,
                       mat_dir_prefix,
                       coupon_dir_prefix,
                       debt_at_par_cp_fn_prefix,
                       eq_fd_cp_fn_prefix)

    dfc = BatchDFColumns(main_params,
                         k_struct_params,
                         fixed_params,
                         debt_vars,
                         equity_vars,
                         share_values,
                         dfcols)


    model_inputs = ModelInputs(Dict{Symbol, Float64}(),
                               ".", ".", ".", ".",
                               ".", ".", ".")


    bp = BatchParams(param_values_dict,
                     common_params,
                     params_order,
                     repeat([[NaN]], inner=ones(Int64, 9)),
                     DataFrame())

    bt = BatchStruct(model,
                     dfn,
                     dfc,
                     model_inputs,
                     bp,
                     coupon_grid,
                     pvb_tol_vec)

    # Set Parameter Combinations
    bt = set_params_combs(bt)

    # Set Directories Paths
    bt = set_main_dir_path(bt)


    return bt
end


# * Batch Jobs - Set & Get Methods
# ** Batch Set Get Params
# include("_batch_set_get_methods/_batch_set_get_params.jl")
# ############## Set Parameter Combinations ##############
function set_param_comb_df(bt)
    # Create Parameter Combinations DataFrame
    pcdf = DataFrame(hcat(bt.bp._params_combs...)')

    cols_dict = Dict([names(pcdf)[i] => Symbol(bt.bp._params_order[i])
                      for i in 1:size(names(pcdf), 1)])

    rename!(pcdf, cols_dict)

    pcdf[!, :comb_num] .= 1:size(pcdf, 1)

    # Important: sort by :m
    sort!(pcdf, :m)
  # Count rows by :m
  # pcdf[!, :m_comb_num] .= by(pcdf, :m, m_comb_num = :m => x -> 1:size(x, 1))[:, :m_comb_num]
  pcdf[!, :m_comb_num] .= combine(groupby(pcdf, ["m"]), :m => (x -> 1:size(x, 1)) => :m_comb_num)[:, :m_comb_num]

    idcols = [:comb_num, :m, :m_comb_num]
    bt.bp.df = pcdf[:, vcat(idcols, [Symbol(x) for x in names(pcdf) if !(Symbol(x) in idcols)])]

    return bt
end


function set_params_combs(bt)
    bt.bp._params_combs = [Array([x1, x2, x3, x4, x5, x6, x7, x8, x9])
                       for x1 in bt.bp._param_values_dict[bt.bp._params_order[1]],
                           x2 in bt.bp._param_values_dict[bt.bp._params_order[2]],
                           x3 in bt.bp._param_values_dict[bt.bp._params_order[3]],
                           x4 in bt.bp._param_values_dict[bt.bp._params_order[4]],
                           x5 in bt.bp._param_values_dict[bt.bp._params_order[5]],
                           x6 in bt.bp._param_values_dict[bt.bp._params_order[6]],
                           x7 in bt.bp._param_values_dict[bt.bp._params_order[7]],
                           x8 in bt.bp._param_values_dict[bt.bp._params_order[8]],
                           x9 in bt.bp._param_values_dict[bt.bp._params_order[9]]]

    return set_param_comb_df(bt)
end
# ########################################################


# ################### Store Parameters ###################
# function set_par_dict(bt, comb_num)
#     bt.mi._svm_dict = Dict()
#     for par in bt._params_order
#         bt.mi._svm_dict[par] = bt._params_combs[comb_num][findall(x-> x==par, bt._params_order)[1]]
#     end

#     bt.mi._svm_dict = merge(bt._common_params, bt.mi._svm_dict)

#     return bt
# end


function set_par_dict(bt; comb_num::Integer=0,
                      m::Float64=NaN, m_comb_num::Integer=0,
                      display_msgs=true)
    # No combination entered
    cond1 = .&(comb_num == 0, (isnan(m) || m_comb_num == 0))
    # Two combinations entered
    cond2 = .&(comb_num > 0, !isnan(m),
               m_comb_num > 0)
    if cond1
        println("Please enter a combination number. Exiting...")
        return bt
    elseif cond2
        println("Two identifiers entered: (i) unique combination ID,
                        (ii) (m, m_comb_num) ID pair.")
        println("Function will use the unique combination ID.")
    end

    bt.mi._svm_dict = Dict()
    if comb_num > 0
        if display_msgs
            println("Setting parameter dictionary using unique combination ID...")
        end
        for par in bt.bp._params_order
            bt.mi._svm_dict[par] = bt.bp._params_combs[comb_num][
                    findall(x-> x==par, bt.bp._params_order)[1]]
        end
    elseif .&(m > 0., m_comb_num > 0)
        if display_msgs
            println("Setting parameter dictionary using (m, m_comb_num) ID pair...")
        end
        # Slice Parameter Combinations DataFrame
        cond = .&(abs.(bt.bp.df[:, :m] .- m) .< 1e-6,
                  abs.(bt.bp.df[:, :m_comb_num] .- m_comb_num) .< 1e-6)
        for par in bt.bp._params_order
            bt.mi._svm_dict[par] = bt.bp.df[cond, par][1]
        end
    end

    bt.mi._svm_dict = merge(bt.bp._common_params, bt.mi._svm_dict)

    return bt
end


function get_par_dict(bt)
    return bt.mi._svm_dict
end


# ########################################################


function get_batch_comb_num(bt; display_msgs::Bool=true,
                            svm_dict::Dict{Symbol,Float64}=Dict{Symbol,Float64}(),
                            iota::Float64=NaN,
                            kappa::Float64=NaN,
                            lambda::Float64=NaN,
                            sigmah::Float64=NaN)

    # Check if all parameters are NaN ########
    nancond = .&(isnan(iota), isnan(kappa), isnan(lambda), isnan(sigmah))
    if .&(isempty(svm_dict), isempty(bt.mi._svm_dict))
        println("Parameter dictionary not found. Returning...")
        return
    # elseif .&(isempty(svm_dict), nancond)
    #     println("Please define at least one parameter value. Returning...")
    #     return
    end
    # ########################################

    # Set Parameter Values: ##################
    if isempty(svm_dict)
        if display_msgs
            println("Setting parameter dictionary to batch object's parameter dictionary")
        end
        svm_dict = bt.mi._svm_dict
    end

    if !nancond
        if display_msgs
            println("Setting parameter values... ")
        end
        if !isnan(iota)
            svm_dict[:iota] = eval(iota)
        end
        if !isnan(kappa)
            svm_dict[:kappa] = eval(kappa)
        end
        if !isnan(lambda)
            svm_dict[:lambda] = eval(lambda)
        end
        if !isnan(sigmah)
            svm_dict[:sigmah] = eval(sigmah)
        end
    end
    # ########################################

    # Find row in the Parameter Combinations DataFrame
    svm_cols = [:lambda, :sigmah]
    cols = [x for x in Symbol.(names(bt.bp.df)) if !(x in vcat([:comb_num, :m_comb_num], svm_cols))]
    LL = [abs.(bt.bp.df[:, col] .- svm_dict[col]) .< 1e-5 for col in cols]
    for col in svm_cols
        if isnan(svm_dict[col])
            LL = vcat(LL, [isnan.(bt.bp.df[:, col])])
        else
            LL = vcat(LL, [abs.(bt.bp.df[:, col] .- svm_dict[col]) .< 1e-5])
        end
    end

    return DataFrame(bt.bp.df[.&(LL...), [:comb_num, :m, :m_comb_num]])
end


function get_batch_comb_numbers(bt, pardict::Dict{Symbol,Array{Float64,1}};
                                idcols::Array{Symbol,1}=[:comb_num, :m, :m_comb_num])
    cond = .&([any.([(abs.(x .- pardict[var]).<1e-6) for x in bt.bp.df[:, var]])
            for var in keys(pardict)]...)


   return bt.bp.df[cond, unique(vcat(idcols, [k for k in keys(pardict)]))]
end


# ** Batch Set Get Dirs
# include("_batch_set_get_methods/_batch_set_get_dirs.jl")
function str_format_fun(a::String,b::Float64)
    return @eval @sprintf($a, $b)
end


# ############## Set Directories ##############
function form_main_dir_path(main_dir::String)
    return joinpath(pwd()[1:findfirst("artur", pwd())[end]],
                    "BondPricing", main_dir)
end


function set_main_dir_path(bt)
    # Main directory is called "Julia". Find it!
    # bt.mi.main_dir_path = joinpath(pwd()[1:findfirst("artur", pwd())[end]],
    #                                "BondPricing", bt.dfn.main_dir)
    bt.mi.main_dir_path = form_main_dir_path(bt.dfn.main_dir)
    return bt
end


function set_batch_res_dir(bt)
    # Path to the Results Directory
    cond = .&(bt.model == "svm",
              abs(bt.mi._svm_dict[:sigmah] - bt.mi._svm_dict[:sigmal]) < 1e-4)
    if cond
        bt.mi.batch_res_dir = string(bt.dfn.res_dir, "Tests")
    else
        bt.mi.batch_res_dir = string(bt.dfn.res_dir, "/", uppercase(bt.model))
    end

    return bt
end


function set_maturity_dir(bt)
    bt.mi.maturity_dir =string(bt.dfn.mat_dir_prefix,
                               str_format_fun(comb_folder_dict[:m][2],
                                              get_par_dict(bt)[:m]))
    return bt
end


function set_comb_res_dir(bt)
    folder_name = string(comb_folder_dict[:mu_b][1], str_format_fun(comb_folder_dict[:mu_b][2],
                                                                    get_par_dict(bt)[:mu_b]),
                         comb_folder_dict[:xi][1], str_format_fun(comb_folder_dict[:xi][2],
                                                                  get_par_dict(bt)[:xi]),
                         comb_folder_dict[:kappa][1], str_format_fun(comb_folder_dict[:kappa][2],
                                                               1e4 * get_par_dict(bt)[:kappa]),
                         comb_folder_dict[:kappa][3],
                         comb_folder_dict[:gross_delta][1], str_format_fun(comb_folder_dict[:gross_delta][2],
                                                                     1e4 * get_par_dict(bt)[:gross_delta]),
                         comb_folder_dict[:gross_delta][3],
                         comb_folder_dict[:iota][1], str_format_fun(comb_folder_dict[:iota][2],
                                                              1e4 * get_par_dict(bt)[:iota]),
                         comb_folder_dict[:iota][3])

    # if "sigma" in bt.bp._params_order
    if bt.model == "cvm"
        bt.mi.comb_res_dir = string(folder_name, comb_folder_dict[:sigmal][1],
                                     str_format_fun(comb_folder_dict[:sigmal][2], get_par_dict(bt)[:sigmal]))
    else
        bt.mi.comb_res_dir = string(folder_name,
                                    comb_folder_dict[:lambda][1], str_format_fun(comb_folder_dict[:lambda][2],
                                                                           get_par_dict(bt)[:lambda]),
                                    comb_folder_dict[:sigmal][1], str_format_fun(comb_folder_dict[:sigmal][2],
                                                                           get_par_dict(bt)[:sigmal]),
                                    comb_folder_dict[:sigmah][1], str_format_fun(comb_folder_dict[:sigmah][2],
                                                                           get_par_dict(bt)[:sigmah]))
    end

    return bt
end


# ############## Get Directories ##############
function get_batch_res_dir(bt)
    return bt.mi.batch_res_dir
end

function get_maturity_dir(bt)
    return bt.mi.maturity_dir
end

function get_comb_res_dir(bt)
    return bt.mi.comb_res_dir
end

# ############## Set Paths ##############
function get_main_dir_path(bt)
    return bt.mi.main_dir_path
end


function set_comb_res_paths(bt)

    # Path to Main Directory
    bt = set_main_dir_path(bt)

    # Generate Folder Names:
    bt = set_batch_res_dir(bt)
    bt = set_maturity_dir(bt)
    bt = set_comb_res_dir(bt)

    # ############ Set Paths ###########
    # Main Results
    bt.mi.batch_res_path = joinpath(get_main_dir_path(bt), get_batch_res_dir(bt))

    # Maturity
    bt.mi.maturity_path = joinpath(bt.mi.batch_res_path, get_maturity_dir(bt))

    # Parameter Combination
    bt.mi.comb_res_path = joinpath(bt.mi.maturity_path, get_comb_res_dir(bt))

    return bt
end


function get_coupon_res_path(bt, coupon::Float64)
    bt = set_comb_res_paths(bt)

    return joinpath(bt.mi.comb_res_path,
                    string(bt.coupon_dir_prefix, str_format_fun(comb_folder_dict[:c][2], coupon)))
end



# ############## Create Directories ##############
function mk_comb_res_dirs(bt)
    # Generate Paths
    bt = set_comb_res_paths(bt)

    # Main Results
    if occursin("Tests", bt.mi.batch_res_path)
        pos = findfirst("Tests", bt.mi.batch_res_path)[1] - 1
        path = bt.mi.batch_res_path[1:pos]
        if !isdir(path)
            mkdir(path)
        end
    end

    main_res_path = split(bt.mi.batch_res_path, uppercase(bt.model))[1]
    if !isdir(main_res_path)
        mkdir(main_res_path)
    end

    if !isdir(bt.mi.batch_res_path)
        mkdir(bt.mi.batch_res_path)
    end

    # Maturity
    if !isdir(bt.mi.maturity_path)
        mkdir(bt.mi.maturity_path)
    end

    # Parameter Combination
    if !isdir(bt.mi.comb_res_path)
        mkdir(bt.mi.comb_res_path)
    end

    return bt
end

function mk_coupon_dir(bt, coupon::Float64)
    mk_comb_res_dirs(bt)
    coupon_res_dir = get_coupon_res_path(bt, coupon)

    # Create coupon sub-folder:
    if !isdir(coupon_res_dir)
        mkdir(coupon_res_dir)
    end

    return coupon_res_dir
end


# ############## Get Directories ##############
function get_results_path(bt, m::Float64)
    set_main_dir_path(bt)
    if bt.mi._svm_dict == nothing
        bt.mi.batch_res_dir = bt.dfn.res_dir
    else
        set_batch_res_dir(bt)
    end

    return string(bt.mi.main_dir_path, "/",
                  bt.mi.batch_res_dir,
                  string(bt.dfn.mat_dir_prefix, str_format_fun(comb_folder_dict[:m][2], m)))
end


# ** Batch Get Objects
# include("_batch_set_get_methods/_batch_get_objs.jl")
# ############## Get Batch & SVM Objects ##############
function get_bt(; model::String="svm", comb_num::Integer=0,
                m::Float64=0., m_comb_num::Integer=0,
                bt=nothing, display_msgs=true)
    # Create Model Object'
    if bt==nothing
        bt = BatchObj(;model=model)
    end

    # Set Combination
    bt = set_par_dict(bt; comb_num=comb_num,
                      m=m, m_comb_num=m_comb_num,
                      display_msgs=display_msgs)

    # Set Directories & Paths (Main, Batch, Maturity, Combination)
    bt = set_comb_res_paths(bt)

    # Set Directories & Paths (Main, Batch, Maturity, Combination)
    bt = set_comb_res_paths(bt)

    return bt
end


function get_bt_cvm(;comb_num::Integer=0,
                    m::Float64=0., m_comb_num::Integer=0,
                    bt=nothing, display_msgs::Bool=true)

    bt = get_bt(;model="cvm",
                comb_num=comb_num,
                m=m, m_comb_num=m_comb_num,
                bt=bt, display_msgs=display_msgs)
    cvm = firm_constructor(bt.mi._svm_dict; model="cvm")

    return bt, cvm
end


function get_bt_svm(;comb_num::Integer=0,
                    m::Float64=0., m_comb_num::Integer=0,
                    bt=nothing, display_msgs::Bool=true,
                    compute_surfs::Bool=false)
    bt = get_bt(;model="svm",
                comb_num=comb_num,
                m=m, m_comb_num=m_comb_num,
                bt=bt, display_msgs=display_msgs)

    # Create Directories
    mk_comb_res_dirs(bt)

    local batch_obj_exists = false
    try
        # batch_obj_exists = check_batch_results(bt)
        # batch_obj_exists = check_batch_file(bt)
        println("Loading SVM object...")
        svm = load_bpr_surfs(bt)

        # Check if Bond Pricing Surfaces exist
        surfs_exist = check_bpr_surfaces(svm)
        batch_obj_exists = true
    catch
        # println("Unable to load Batch object file. Recomputing... ")
        println("Unable to locate Batch object file. Recomputing...")
        batch_obj_exists = false
    end


    println(string("Batch object exists: ", batch_obj_exists))
    # Check if results exist
    if !batch_obj_exists
    #     # Load SVM object
    #     println("Loading SVM object...")
    #     # svm = load_batch_results(bt)["svm"]
    #     svm = load_bpr_surfs(bt)

    #     # Check if Bond Pricing Surfaces exist
    #     surfs_exist = check_bpr_surfaces(svm)
    # else
        # Create SVM object
        svm = firm_constructor(bt.mi._svm_dict; model="svm")
        surfs_exist = false
    end

    if !surfs_exist | compute_surfs
        # Compute Bond Price Inputs
        svm = @time bpr_surfs(svm)

        # Save results
        save_svm_surfaces(svm, bt.mi.comb_res_path)
    end

    # Interpolate Bond Pricing Surfaces
    println("Interpolating bond pricing surfaces...")
    svm = @time bpr_interp(svm)

    return bt, svm
end


function get_bt_mobj(;model::String="svm",
                     comb_num::Integer=0,
                     m::Float64=0., m_comb_num::Integer=0,
                     bt=nothing,
                     display_msgs::Bool=true,
                     compute_surfs::Bool=false)

    if model == "svm"
        return get_bt_svm(;comb_num=comb_num,
                          m=m, m_comb_num=m_comb_num,
                          bt=nothing, display_msgs=display_msgs,
                          compute_surfs=compute_surfs)
    elseif model == "cvm"
        return get_bt_cvm(;comb_num=comb_num,
                          m=m, m_comb_num=m_comb_num,
                          bt=nothing, display_msgs=display_msgs)
    else
        println("Please choose model cvm or svm. Exiting...")
        return
    end
end


# * Batch Jobs - File Methods
# ** Batch Object Files
# include("_batch_file_methods/_batch_obj_files.jl")
function save_svm_surfaces(svm, file_path)
    bond_path_fname = string(file_path, "/bond_surfaces", ".jld")
    bs = Dict("vtgrid" => svm.bs.vtgrid,
              "ttmgrid" => svm.bs.ttmgrid,
              "vbhlgrid" => svm.bs.vbhlgrid,
              "f11_surf" => svm.bs.f11_surf,
              "f12_surf" => svm.bs.f12_surf,
              "f13_surf" => svm.bs.f13_surf,
              "f21_surf" => svm.bs.f21_surf,
              "f22_surf" => svm.bs.f22_surf)
    # save(bond_path_fname, "bs", bs)
    jldopen(bond_path_fname, "w") do file
        write(file, "bs", bs)
    end
end


function save_batch_results(bt; svm=nothing, fname::String=batch_file_name)
    bt = mk_comb_res_dirs(bt)
    path_fname = string(bt.mi.comb_res_path, "/", fname, ".jld")


    if !(svm==nothing)
#        save(path_fname, "bt", bt, "svm", svm)
        jldopen(path_fname, "w") do file
            write(file, "bt", bt)
            write(file, "svm", svm)
        end

        save_svm_surfaces(svm, bt.mi.comb_res_path)
    else
        # save(path_fname, "bt", bt)
        jldopen(path_fname, "w") do file
            write(file, "bt", bt)
        end
    end
end


function check_batch_file(bt; fname::String=batch_file_name)
    path_fname = string(bt.mi.comb_res_path, "/", fname, ".jld")

    return .&(isfile(path_fname), stat(path_fname).size > 0)
end


function load_batch_results(bt; fname::String=batch_file_name)
    path_fname = string(bt.mi.comb_res_path, "/", fname, ".jld")

    exist_cond = false
    if check_batch_file(bt; fname=batch_file_name)
        try
            res = load(path_fname)
            # Check if both bt and svm objects are present:
            exist_cond = .&([x in keys(res) for x in ["bt", "svm"]]...)

            if !exist_cond
                println("Batch or SVM object missing. ")
            end
        catch
            println("Unable to load batch object file.")
        end
    else
        println("Batch object file not found.")
    end

    if exist_cond
        return load(path_fname)
    end
end


function check_batch_results(bt; fname::String=batch_file_name)
    # if batch object exists, check that surfaces have
    # been computed.
    svm_exists = false

    if check_batch_file(bt; fname=fname)
        # Load Batch Object:
        res = load_batch_results(bt)

        # Check that SVM object is in results object:
        return ("svm" in keys(res))
    end

    return svm_exists
end


function check_bpr_surfaces(svm)
    # Check that SVM object contains bond pricing surfaces:
    return !(any(isnan, svm.bs.f11_surf) ||
             any(isnan, svm.bs.f12_surf) ||
             any(isnan, svm.bs.f13_surf) ||
             any(isnan, svm.bs.f21_surf) ||
             any(isnan, svm.bs.f22_surf))
end


function load_bpr_surfs(bt; batch_file_name::String=batch_file_name, interp_surfs::Bool=true)
    try check_batch_results(bt; fname=batch_file_name)
        svm = load_batch_results(bt, fname=batch_file_name)["svm"]
    catch
        println("Batch Object is incompatible! Loading surfaces only instead.")
        svm = firm_constructor(bt.mi._svm_dict)
        svmbs = load(string(bt.mi.comb_res_path, "/bond_surfaces.jld"))["bs"]
        svm.bs.vtgrid = svmbs["vtgrid"]
        svm.bs.ttmgrid = svmbs["ttmgrid"]
        svm.bs.vbhlgrid = svmbs["vbhlgrid"]
        svm.bs.f11_surf = svmbs["f11_surf"]
        svm.bs.f12_surf = svmbs["f12_surf"]
        svm.bs.f13_surf = svmbs["f13_surf"]
        svm.bs.f21_surf = svmbs["f21_surf"]
        svm.bs.f22_surf = svmbs["f22_surf"]
    end

    # Interpolate Bond Pricing Surfaces
    if interp_surfs
        println("Interpolating Bond Pricing Surfaces...")
        svm = @time bpr_interp(svm)
    end

    return svm
end


# ** Batch Equity Files
# include("_batch_file_methods/_batch_equity_files.jl")
function conditionally_load_df(df, path_fname::String,
                               types::Array{DataType,1}=vcat(fill(Float64, 8), [Int64, Int64]))
    try
        return eval(df)
    catch
        return CSV.read(path_fname, DataFrame; types=types)
    end
end


function save_eq_results(bt, res, df::DataFrame, eq_fd_cp_path_fname::String)
    sol = sort!(vcat(df[:, bt.dfc.dfcols], res[:, bt.dfc.dfcols]), [:p])

    # Identify solution
    sol[!, :opt] .= .&(abs.(sol[:, :abs_debt_diff] .-
                         minimum(sol[:, :abs_debt_diff])) .< 1e-4,
                        abs.(sol[:, :eq_deriv]) .< 1e-4)

    CSV.write(eq_fd_cp_path_fname, sol)
end


function load_eq_results(bt, svm, dfn; use_all_eqdf::Bool=true)
    try
      #= coltypes = vcat(fill(Float64, 16), [Bool], fill(Float64, 16)) =#
      if .&(use_all_eqdf,  string("all_", dfn) in readdir(bt.mi.comb_res_path))
        all_eqdf = CSV.read(string(bt.mi.comb_res_path, "/all_", dfn), DataFrame)

        #= if size(all_eqdf, 2) == size(coltypes, 1) =#
        #=   all_eqdf = CSV.read(string(bt.mi.comb_res_path, "/all_", dfn), =#
        #=                       types=coltypes) =#
        #= end =#

       # Compute candidate vb for each p value
        eqdf_final = Batch.eq_fd_processing(bt, svm, all_eqdf)

        # Compute (p, vb) solution for c value
        res = Batch.eq_fd_sol(bt, svm, eqdf_final)

        # Save results
        # Batch.save_eq_results(bt, res, eqdf_final,
        #                      string(bt.mi.comb_res_path, "/", dfn))

        try
          sol = sort!(vcat(eqdf_final[:, bt.dfc.dfcols], res[:, bt.dfc.dfcols]), [:p])
        catch 
          common_cols = [x for x in names(eqdf_final) if (x in names(res))]
          sol = sort!(vcat(eqdf_final[:, common_cols], res[:, common_cols]), [:p])
        end

        # Identify solution
        sol[!, :opt] .= .&(abs.(sol[:, :abs_debt_diff] .-
                            minimum(sol[:, :abs_debt_diff])) .< 1e-4,
                       abs.(sol[:, :eq_deriv]) .< 1e-4)

        # Drop duplicates
        unique!(sol, [:p])

        CSV.write(string(bt.mi.comb_res_path, "/", dfn), sol)

        return sol
      elseif !use_all_eqdf
        sol = CSV.read(string(bt.mi.comb_res_path, "/", dfn), DataFrame)#;
        #               types=vcat(coltypes, [Bool]))

        return unique!(sol, [:p, :vb])
      else
        return DataFrame()
      end
    catch
      return DataFrame()
    end
end


# ** Batch Results Files
# include("_batch_file_methods/_batch_results_files.jl")
function load_svm_opt_results(bt;
                              firm_obj_fun::Symbol=:firm_value,
                              comb_num::Integer=0,
                              m::Float64=NaN,
                              m_comb_num::Integer=0,
                              display_msgs::Bool=true,
                              toldf::DataFrame=Batch.toldf,
                              drop_fail::Bool=false,
                              save_soldf::Bool=true,
                              soldf_name::String=soldf_name,
                              interp_polyOrder::Integer=3,
                              filter_windowSize::Integer=5 * 10^2 + 1,
                              filter_polyOrder::Integer=3,
                              save_optdf::Bool=true,
                              optdf_name::String=optdf_name)

    bt = get_bt(;comb_num=comb_num,
                m=m, m_comb_num=m_comb_num,
                display_msgs=display_msgs)

    opt_k_struct_found = false
    if string(optdf_name, ".csv") in readdir(bt.mi.comb_res_path)
        cols = opt_k_struct_cols(bt)
        # coltypes = map(x -> begin
        #                       if x == :eq_negative
        #                         return Bool
        #                       end
        #                       return Float64
        #                     end, cols)
        optDF = CSV.read(string(bt.mi.comb_res_path, "/", optdf_name, ".csv"), DataFrame;
                         types=opt_k_struct_df_coltypes)

        if !isnan(optDF[1, :c])
            opt_k_struct_found=true
        end
    end

    # if opt_k_struct_found
    #     #optDF = hcat(get_batch_comb_num(bt; display_msgs=display_msgs)[[:comb_num, :m_comb_num]], optDF)
    #     combdf = get_batch_comb_num(bt)
    #     optDF[:comb_num] =combdf[1, :comb_num]
    #     optDF[:m_comb_num] =combdf[1, :comb_num]
    #     id_cols = [:comb_num, :m, :m_comb_num]
    #     cols_order = vcat(id_cols, [x for x in names(optDF) if !(x in id_cols)])
    #     optDF = optDF[cols_order]
    if !opt_k_struct_found
        optDF = get_opt_results(bt;
                                comb_num=comb_num,
                                firm_obj_fun=firm_obj_fun,
                                display_msgs=display_msgs,
                                toldf=toldf,
                                drop_fail=drop_fail,
                                save_soldf=save_soldf,
                                soldf_name=soldf_name,
                                interp_polyOrder=interp_polyOrder,
                                filter_windowSize=filter_windowSize,
                                filter_polyOrder=filter_polyOrder,
                                save_optdf=save_optdf,
                                optdf_name=optdf_name)
    end
    return optDF
end


function load_svm_opt_results_df(bt;
                                 firm_obj_fun::Symbol=:firm_value,
                                 m::Float64=NaN,
                                 opt_k_struct_df_name::String=opt_k_struct_df_name,
                                 coltypes::Array{DataType,1}=opt_k_struct_df_coltypes,
                                 recompute::Bool=false)
    if !isnan(m)
        bt = get_bt(;bt=bt, m=m, m_comb_num=1)
    end

    extension = lowercase(string(firm_obj_fun))
    dfname = string(opt_k_struct_df_name, "_", extension)

    if recompute || !(string(dfname, ".csv") in readdir(bt.mi.maturity_path))
        println("Compiling results...")
        optDF = compile_svm_opt_results(bt; m=m,
                                        firm_obj_fun=firm_obj_fun,
                                        opt_k_struct_df_name=opt_k_struct_df_name)

    else
        if isnan(m)
            fpath = bt.mi.batch_res_path
        else
            fpath = bt.mi.maturity_path
        end
        println("Loading optimal results dataframe...")
        optDF = CSV.read(string(fpath, "/", dfname, ".csv"), DataFrame; types=coltypes)
    end

    # if !(:MBR in names(optDF))
    #     optDF[:MBR] = get_mbr(unique(optDF[:V0])[1], optDF[:debt], optDF[:equity])
    # end
    return optDF
end



# * Batch Jobs - Debt at Par and Equity Finite Differences Method
# ** Batch PVB Grids
# include("_batch_core/_batch_pvb_grids.jl")
function get_cvm_p_debt_diff(svm, sigma, mu_b, c, pgrid; N=10000)
    # Compute CVM VBs
    cvm_vb = [get_cvm_vb(svm, sigma; mu_b=mu_b, c=c, p=p) for p in pgrid]

    # Interpolate CVM VBs
    cvm_vb_fun = Dierckx.Spline1D(pgrid, cvm_vb, k=3, bc="extrapolate")

    # Compute Debt Values
    cvm_debt = fetch(@spawn [get_cvm_debt_price(svm, cvm_vb_fun(p), sigma;
                                            mu_b=mu_b, c=c, p=p) for p in pgrid])

    # Compute Debt-P
    aggPgrid = [get_agg_p(svm; mu_b=mu_b, p=p) for p in pgrid]
    cvm_debt_diff = cvm_debt - aggPgrid

    # Interpolate and Find Optimal Value
    pgrid_ref = range(pgrid[1], stop=pgrid[end], length=N)
    cvm_dff = Dierckx.Spline1D(pgrid, cvm_debt_diff, k=3, bc="extrapolate")
    # Get highest pj8 value:
    cvm_pOpt = reverse(pgrid_ref)[argmin(abs.(reverse(cvm_dff(pgrid_ref))))]
    cvm_vbOpt = cvm_vb_fun(cvm_pOpt)

    return Dict("cvm_vb_fun" => cvm_vb_fun,
                "cvm_debt" => cvm_debt,
                "cvm_dff" => cvm_dff,
                "cvm_pOpt" => cvm_pOpt,
                "cvm_vbOpt" => cvm_vbOpt)
end


function p_candidates(svm; pmin=5, pmax=NaN, mu_b=NaN, c=NaN, N1=20, N2=10000, lower_adj=.85, upper_adj=1.15)
    if isnan(pmax)
        pmax = .9 * svm.pm.V0
    end

    if isnan(mu_b)
        mu_b = svm.mu_b
    end

    if isnan(c)
        c = svm.c
    end

    # Form p grids
    pgrid = range(pmin, stop=pmax, length=N1)

    # Get VB, Debt, Debt Diff, Debt at Par values
    cvml_dict = get_cvm_p_debt_diff(svm, svm.pm.sigmal, mu_b, c, pgrid; N=N2)
    cvmh_dict = get_cvm_p_debt_diff(svm, svm.pm.sigmah, mu_b, c, pgrid; N=N2)

    # Form Adjusted pgrid
    p_lower = lower_adj * minimum([cvml_dict["cvm_pOpt"], cvmh_dict["cvm_pOpt"]])
    p_upper = upper_adj * maximum([cvml_dict["cvm_pOpt"], cvmh_dict["cvm_pOpt"]])
    adj_pgrid = range(p_lower, stop=p_upper, length=N1)

    # Ensure VB values are within bounds

    return Dict("cvml" => cvml_dict, "cvmh" => cvmh_dict, "adj_pgrid" => adj_pgrid)
end


function bounded_vbs(svm, cont, p)

    cvm_vbl = cont["cvml"]["cvm_vb_fun"](p)
    cvm_vbh = cont["cvmh"]["cvm_vb_fun"](p)

    # Form VB grid
    vblmin = .5 * minimum([cvm_vbl, cvm_vbh])
    vblmax = 1.5 * maximum([cvm_vbl, cvm_vbh])
    vblgrid = range(vblmin, stop=vblmax, length=20)

    # Compute vbh/vbl ratios
    ratio = [cvm_vbh/vbl for vbl in vblgrid]

    # Interpolate ratios
    ratiof = Dierckx.Spline1D(vblgrid, ratio, k=3, bc="extrapolate")

    # Refine vblgrid
    vblgrid_ref = range(vblmin, stop=vblmax, length=10000)

    # Ensure vbl values satisfy bhl_min <= cvm_vbh/vbl <= bhl_max
    bounded_vblgrid = [vb for vb in vblgrid_ref if
                        (ratiof(vb) <= svm.bi.vbhlmax) & (ratiof(vb) >= svm.bi.vbhlmin)]
    vblmin = minimum(bounded_vblgrid)
    vblmax = maximum(bounded_vblgrid)

    return vblmin, vblmax
end


# ** Batch Debt at Par Functions
# include("_batch_core/_batch_debt_at_par_funs.jl")
function get_cvm_vb_funs(svm, c, p_grid)
    cvm_vbl = [AnalyticFunctions.get_cvm_vb(svm, svm.pm.sigmal; c=c, p=p) for p in p_grid]
    cvm_vbh = [AnalyticFunctions.get_cvm_vb(svm, svm.pm.sigmah; c=c, p=p) for p in p_grid]
    cvm_vbl_fun = Dierckx.Spline1D(p_grid, cvm_vbl, k=3, bc="extrapolate")
    cvm_vbh_fun = Dierckx.Spline1D(p_grid, cvm_vbh, k=3, bc="extrapolate")

    return cvm_vbl_fun, cvm_vbh_fun
end

# #######################################################
# ##################### DEBT AT PAR #####################
# #######################################################
function debt_at_par_grid_bounds(svm, c, x_min, x_max; p=NaN, vbl=NaN,
                                 f_low=.75, f_high=1.25)
    xgrid = range(f_low * x_min, stop=f_high * x_max, length=50)


    vbh = get_cvm_vb(svm, svm.pm.sigmah; c=c, p=p)
    if isnan(vbl)
        vbhl_vec = [vbh/xval for xval in xgrid]
    elseif isnan(p)
        vbhl_vec = [get_cvm_vb(svm, svm.pm.sigmah;
                               c=c, p=x)/vbl for x in xgrid]
    else
        println("Please enter a value for p or vbl")
        return
    end

    # Interpolate
    vbhlf = Dierckx.Spline1D(xgrid, vbhl_vec, k=3, bc="extrapolate")

    # Refine
    xgrid_ref  = range(f_low * x_min, stop=f_high * x_max, length=10^5)
    vbhl_ref = [vbhlf(xval) for xval in xgrid_ref]

    # Bounds
    x0 = xgrid_ref[argmin([abs(vbhl - 1.025 * svm.bi.vbhlmin) for vbhl in vbhl_ref])]
    x1 = xgrid_ref[argmin([abs(vbhl - .975 * svm.bi.vbhlmax) for vbhl in vbhl_ref])]


    # Now Make Sure vt bounds are satisfied:
    # xmax = minimum([maximum([x0,x1]),  .975 * svm.pm.V0 * exp(- minimum(svm.bs.vtgrid))])

    return minimum([x0,x1, .95 * svm.pm.V0]), minimum([maximum([x0,x1]), svm.pm.V0]) # 1.001 * minimum([x0, x1]), .999 * maximum([x0, x1])
end


function debt_at_par_check_conds(svm, df, xvar, xmin, xmax; disp_option=false)
    # DF should satisfy the following conditions:
    # 1. At least 5 (P, VB) pairs with precision <= 2*1e-2 (2% off-mark)
    # 2. At least 2 out of 5 (P, VB) pairs should be in the interval:
    #                       [.9 * p_min, 1.1 * p_max)]
    # 3. At least 2 out of 5 pairs should have precision <= 1e-2
    cond1 = sum(df[:, :abs_debt_per_diff] .<= 2e-2) >= 5
    cond2 = sum(.&(df[:, xvar] .>= xmin * .9, df[:, xvar] .<= xmax * 1.1)) >= 2
    cond3 = sum(df[:, :abs_debt_per_diff] .<= 1e-2) >= 2
    conds = cond1 & cond2 & cond3

    if disp_option
        println(string("Condition 1: ", cond1))
        println(string("Condition 2: ", cond2))
        println(string("Condition 3: ", cond3))
        println(string("Conds: ", conds))
    end

    return conds
end


function debt_at_par_diffs(svm, df, xvar; p=NaN, sort_vars=[])
    if xvar ==:p
        aggP = [get_agg_p(svm, p=p) for p in df[:, :p]]
    elseif !isnan(p)
        aggP = get_agg_p(svm, p=p)
    else
        println("Missing p value(s). Exiting...")
    end

    df[!, :debt_diff] .= df[:, :debt] .- aggP
    df[!, :abs_debt_diff] .= abs.(df[:, :debt_diff])
    df[!, :debt_per_diff] .= (df[:, :debt_diff]) ./ aggP
    df[!, :abs_debt_per_diff] .= abs.(df[:, :debt_per_diff])

    if isempty(sort_vars)
        return sort!(df, xvar)
    else
        return sort!(df, sort_vars)
    end
end


function debt_at_par_new_df(svm, c, p, df, xvar, debt_vars=debt_vars)
    xgrid_new = (Array(df[xvar])[2:end] + Array(df[xvar])[1:end-1])./ 2.
    df_new = DataFrame(hcat(reverse(xgrid_new), 999. * (ones(length(xgrid_new), size(df, 2) - 1))),
                       vcat([xvar, :debt], debt_vars))
    df_new[!, :debt] .= fetch(@spawn [get_svm_debt_price(svm, vbl; c=c, p=p)
                                  for vbl in df_new[:, xvar]])

    return sort!(vcat(df, df_new), :vb)
end


function debt_at_par_main(svm, c, df, xvar, xgrid_ref;
                          vbl=NaN, p=NaN, debt_vars=debt_vars, disp_option=true)
    # Notice here I use the percentage difference, not
    # the absolute percentage difference. I want to see
    # where the function crosses the zero line.

    # Interpolate
    if xvar == :p
        agg_p_vec = [get_agg_p(svm, p=p) for p in df[:, :p]]
        diff_interp_fun = Dierckx.Spline1D(df[:, :p],
                              (df[:, :debt] .- agg_p_vec) ./ agg_p_vec,
                                           k=3, bc="extrapolate")
    else
        diff_interp_fun = Dierckx.Spline1D(df[:, :vb],
                              (df[:, :debt] .- get_agg_p(svm, p=p)) ./ get_agg_p(svm, p=p),
                                           k=3, bc="extrapolate")
    end

    diffs = diff_interp_fun(xgrid_ref)

    # Find Candidate for the Global Minimum:
    xroot = reverse(xgrid_ref)[argmin(abs.(reverse(diffs)))]
    if disp_option
        println(string(xvar, " root: ", xroot))
    end

    if sum(abs.(df[:, xvar] .- xroot) .< 1e-4) == 0
        if xvar == :p
            debt_root = get_svm_debt_price(svm, vbl; p=xroot)
        else
            debt_root = get_svm_debt_price(svm, xroot; c=c, p=p)
        end

        sort!(push!(df, [xroot, debt_root, NaN, NaN, NaN, NaN]),
              vcat([xvar, :debt], debt_vars))
    end

    return debt_at_par_diffs(svm, df, xvar; p=p)
end


# ATTENTION: NEED TO ADJUST P TO ACCOUNT FOR DIFFERENT M!!!!
# ATTENTION: p_min = np.min([cvml_p,cvmh_p])
# ATTENTION: p_max = 1.5 * np.max([cvml_p,cvmh_p])
function svm_debt_at_par(svm, x_min, x_max, c; p=NaN, vbl=NaN,
                         N1=15, N2=10^5, f_low=.85, f_high=1.15,
                         debt_vars=debt_vars, disp_option=true)

    start_tic = time_ns()

    # Initial Check:
    if isnan(p) & isnan(vbl)
        println("Please enter a value for p or vbl")
        return
    end

    # Form Grids:
    xvar = [:p, :vb][isnan.([p, vbl])][1]
    fixed_var = [:p, :vb][isnan.([p, vbl]).==false][1]
    fixed_var_value = [p, vbl][isnan.([p, vbl]).==false][1]
    xmin, xmax = debt_at_par_grid_bounds(svm, c, x_min, x_max;
                                         p=p, vbl=vbl, f_low=.75, f_high=1.25)
    xgrid = range(xmin, stop=xmax, length=N1)
    xgrid_ref = range(xmin, stop=xmax, length=N2)

    # Compute Debt Values
    df = DataFrame(hcat(reverse(xgrid), 999. * ones(length(xgrid), 5)),
                        vcat([xvar, :debt], debt_vars))

    if xvar == :p
        df[!, :debt] .= fetch(@spawn [get_svm_debt_price(svm, vbl; c=c, p=pval) for pval in df[:, :p]])
    else
        df[!, :debt] .= fetch(@spawn [get_svm_debt_price(svm, vblval; c=c, p=p) for vblval in df[:, :vb]])
    end

    # Sort DataFrame
    sort!(df, xvar)

    # Compute Candidate Solution:
    df = debt_at_par_main(svm, c, df, xvar, xgrid_ref,
                          vbl=vbl, p=p, disp_option=disp_option)

    # Check Conditions:
    conds = debt_at_par_check_conds(svm, df, xvar, x_min, x_max,
                                    disp_option=disp_option)

    count = 0
    if disp_option
        println(string("Counter: ", count))
        println("Unique ", xvar, " values: ", length(unique(df[:, xvar])))
    end

    while !conds & (length(unique(df[:, xvar])) < 20) & (count < 10)
        if disp_option
            println("While loop")
        end

        # Add Data Points:
        df = debt_at_par_new_df(svm, c, p, df, xvar)

        # Compute Candidate Solution:
        df = debt_at_par_main(svm, c, df, xvar,
                              xgrid_ref, vbl=vbl, p=p)

        # Check Conditions:
        conds = debt_at_par_check_conds(svm, df, xvar, x_min, x_max,
                                        disp_option=disp_option)

        # Update Counter:
        count += 1

        if disp_option
            println(string("Counter: ", count))
            println("Unique ", xvar, " values: ", length(unique(df[xvar])))
        end
    end

    if disp_option
        println("Debt at par: preparing results...")
    end
    df[!, :c] .= c
    df[!, fixed_var] .= fixed_var_value
    df[!, :count] .= count
    df[!, :time] .= time_ns() - start_tic

    cols = vcat([:c, fixed_var, xvar, :debt],
                debt_vars, [:count, :time])
    df = df[:, cols]

    if disp_option
        println("returning results...")
        println(string("length: ", size(df, 1)))
    end

    return df
end
# #######################################################


# ** Batch Equity Finite Differences
# include("_batch_core/_batch_eq_fd.jl")
function df_slicer(df, p)
   return df[abs.((df[:, :p] .- p)) .< 1e-4, :]
end


function abs_debt_diff_cond(df, p; tol=.05)
    # Slice DataFrame
    return minimum(df_slicer(df, p)[:, :abs_debt_diff]) < tol
end


function filter_debt_values(bt, svm, df; N1=50, N2=10)
    pval = unique(df[:, :p])[1]

    debtf = Dierckx.Spline1D(df[:, :vb], df[:, :debt], k=3, bc="extrapolate")
    vbgrid = range(minimum(df[:, :vb]), stop=maximum(df[:, :vb]), length=N1)
    aggP = get_agg_p(svm, p=pval)
    debt_diff = debtf(vbgrid) .- aggP

    # Find VB values for which the absolute debt diff is the smaller:
    pos = partialsortperm(abs.(debt_diff), 1:N2)
    vbvals = vbgrid[pos]

    sdf = DataFrame(vb=vbvals, debt=debtf(vbvals))
    sdf[!, :c] .= unique(df[:, :c])[1]
    sdf[!, :p] .= unique(df[:, :p])[1]
    sdf = debt_at_par_diffs(svm, sdf, :vb, p=pval)

    return sdf[:, vcat([:c, :p, :vb, :debt], bt.dfc.debt_vars)]
end


function filter_batch_I_df(bt, svm, df; tol=.05, N1=50, N2=10)
    LL = fetch(@spawn [filter_debt_values(bt, svm, df_slicer(df, p); N1=N1, N2=N2)
                       for p in unique(df[:, :p])])
    sdf=vcat(LL...)

    # Get Filtered pgrid
    pgrid =  [p for p in unique(sdf[:, :p]) if
              (minimum(abs.(df_slicer(sdf, p)[:, :debt_diff])) < tol)]
    while size(pgrid, 1) < 5
        tol = 1.25 * tol
        pgrid =  [p for p in unique(sdf[:, :p]) if
                  (minimum(abs.(df_slicer(sdf, p)[:, :debt_diff])) < tol)]
    end

    # Return Filtered DataFrame
    return df[findall(in(pgrid), df[:, :p]), :]
end


function eq_fd_method(bt, svm, df; mu_b=NaN, c=NaN)
    if isnan(mu_b)
        mu_b=svm.mu_b
    end

    if isnan(c)
        c=svm.c
    end

    res = @time fetch(@spawn [eq_fd(svm; vbl=df[i, :vb], mu_b=mu_b, c=c, p=df[i, :p])
                              for i in 1:size(df, 1)])

    # Collect Results
    tmp = vcat(res...)

    # Compute Debt Differences
    eqdf_all = debt_at_par_diffs(svm, tmp, :p; sort_vars=[:p, :vb])

    # Rearrange Columns
    return eqdf_all[:, bt.dfc.dfcols]
end


function interp_values(res, df, xvar::Symbol, interp_cols::Array{Symbol,1};
                       k::Int64=3, bc::String="extrapolate")
    for col in interp_cols
        ffun = Dierckx.Spline1D(df[:, xvar], df[:, col]; k=k, bc=bc)
        res[!, col] .= ffun(res[:, xvar])
    end
    return res
end


function non_interp_values(svm, df)
    # Debt Values
    df[!, :abs_debt_diff] .= abs.(df[:, :debt_diff])
    df[!, :abs_debt_per_diff] .= abs.(df[:, :debt_per_diff])

    # Equity
    df[!, :eq_negative] .= (df[:, :eq_min_val] .< -.005)

    # Share Values
    df[!, :firm_value] .= df[:, :debt] .+ df[:, :equity]
    df[!, :leverage] .= (df[:, :debt]./df[:, :firm_value]) .* 100
    df[!, :MBR] .= (df[:, :equity]./(svm.pm.V0 .- df[:, :debt]) .- 1) .* 100

    return df
end


function eq_deriv_root_search(svm, df, p; mu_b=NaN, c=NaN, N=10^5)
    if isnan(mu_b)
        mu_b=svm.mu_b
    end

    if isnan(c)
        c=svm.c
    end

    # Create DataFrame
    res = DataFrame(:p => p)

    # Filter DataFrame
    fdf = df[abs.(df[:, :p] .- p) .< 1e-4, :]

    # Interpolate Values
    eq_deriv_fun = Dierckx.Spline1D(fdf[:, :vb], fdf[:, :eq_deriv], k=3, bc="extrapolate")

    # Compute optimal VB:
    vbroots = roots(eq_deriv_fun; maxn=8)
    if !isempty(vbroots)
        # debt_interp = Dierckx.Spline1D(fdf[:vb], fdf[:debt], k=3, bc="extrapolate")
        # abs_debt_diff = abs.(debt_interp(vbroots) .- get_agg_p(svm, p=p))
        # res[:vb] = vbroots[argmin(abs_debt_diff)]
        eq_min_val_interp = Dierckx.Spline1D(fdf[:, :vb], fdf[:, :eq_min_val], k=3, bc="extrapolate")
        abs_eq_min_val = abs.(eq_min_val_interp(vbroots))
        res[!, :vb] .= vbroots[argmin(abs_eq_min_val)]
    else
        ref_vbgrid = range(minimum(fdf[:, :vb]), stop=maximum(fdf[:, :vb]), length=N)
        res[!, :vb] .= ref_vbgrid[argmin(abs.(eq_deriv_fun(ref_vbgrid)))]
    end

    # Equity Values
    res[!, :eq_deriv] .= eq_deriv_fun(res[:, :vb])

    # Interpolate Functions
    interp_cols = vcat([:debt, :equity],
                       [:eq_deriv_min_val,
                        :eq_min_val, :eq_vb])
    res = interp_values(res, fdf, :vb, interp_cols)

    # Debt Values
    res = debt_at_par_diffs(svm, res, :p)

    # Fixed Values
    return non_interp_values(svm, res)

    # # Share Values
    # res[:firm_value] = res[:debt] .+ res[:equity]
    # res[:leverage] = (res[:debt]./res[:firm_value]) .* 100
    # res[:ROE] = (res[:equity]./(svm.pm.V0 .- res[:debt]) .- 1) .* 100.

    # return res
end


function eq_fd_processing(bt, svm, df; mu_b=NaN, c=NaN, N=10^5)
    if isnan(mu_b)
        mu_b=svm.mu_b
    end

    if isnan(c)
        c=svm.c
    end

    params_dict = Dict()
    for par in vcat(bt.dfc.main_params, [:mu_b, :m, :c], bt.dfc.fixed_params)
        params_dict[par] = unique(df[:, par])[1]
    end

    # Compute the Solutions for each pair (c, p)
    tmp = fetch(@spawn [eq_deriv_root_search(svm, df, p; mu_b=mu_b, c=c, N=N)
                        for p in unique(df[:, :p])])

    # Collect results
    eqfinal = hcat(vcat(tmp...), repeat(DataFrame(params_dict), outer=size(tmp, 1)))

    # Rearrange columns
    try
      return eqfinal[:, bt.dfc.dfcols]
    catch
      return eqfinal
    end
end


function eq_fd_sol(bt, svm, df; N=10^5)
    # Find p for which |D - P| = 0:
    ddiff = Dierckx.Spline1D(df[:, :p], df[:, :debt_diff], k=3, bc="extrapolate")
    pgrid = range(minimum(df[:, :p]), stop=maximum(df[:, :p]), length=N)
    res = DataFrame(:p => pgrid[argmin(abs.(ddiff(pgrid)))])
    res[!, :debt_diff] .= ddiff(res[:, :p])

    # Interpolate
    interp_cols = vcat([:vb, :debt, :equity],
                       [:debt_per_diff],
                       [:eq_deriv, :eq_deriv_min_val,
                        :eq_min_val, :eq_vb])
    res = interp_values(res, df, :p, interp_cols)
    res = non_interp_values(svm, res)

    # Get Non-Variable Values
    cols = vcat(bt.dfc.main_params, [:mu_b, :m, :c], bt.dfc.fixed_params)
    for col in cols
        res[!, col] .= df[1, col]
    end

    try 
      return res[:, bt.dfc.dfcols]
    catch
      return res
    end
end


# * Processing Results
# ** Batch Diagnosis
# include("_batch_results/_batch_diagnosis.jl")
function diag_df(bt, i::Integer)
    combdf = DataFrame(bt.bp.df[i, :])

    bt = Batch.set_par_dict(bt; comb_num=combdf[1, :comb_num], display_msgs=false)
    bt = Batch.set_comb_res_paths(bt)

    combdf = hcat(combdf, DataFrame(folder_name = bt.mi.comb_res_dir,
                                    folder_created = isdir(bt.mi.comb_res_path),
                                    batch_obj=false,
                                    modified = " ",
                                    bond_surfs = false,
                                    count = 0))

    # Count solutions:
    if combdf[1, :folder_created]
        dfs_list = [x for x in readdir(bt.mi.comb_res_path) if
                    occursin(bt.dfn.eq_fd_cp_fn_prefix, x)]
        combdf[!, :count] .= size(dfs_list, 1)
        combdf[!, :batch_obj] .= sum([occursin(Batch.batch_file_name, x) for
                                  x in readdir(bt.mi.comb_res_path)]) > 0

        if combdf[1, :batch_obj]
            batch_obj_file = string(bt.mi.comb_res_path, "/",
                                    Batch.batch_file_name, ".jld")
            combdf[!, :modified] .= string(Dates.unix2datetime(stat(batch_obj_file).mtime))
        end
        combdf[!, :bond_surfs] .= sum([occursin("bond_surfaces", x) for
                                   x in readdir(bt.mi.comb_res_path)]) > 0
    end

    return combdf
end


function diagnosis(bt)
    # combDFs = fetch(@spawn [diag_df(bt, comb) for comb in
    #                         range(1, stop=size(bt._params_combs[:], 1))])

    combDFs = fetch(@spawn [diag_df(bt, i) for i in 1:nrow(bt.bp.df)])
    cols1 = [:comb_num, :m, :m_comb_num, :folder_created, :batch_obj,
             :modified, :bond_surfs, :count, :folder_name]
    cols2 = [x for x in Symbol.(names(combDFs[1])) if !(x in cols1)]
    return sort!(vcat(combDFs...)[:, vcat(cols1, cols2)], [:m, :m_comb_num])
end


# ** Batch Process Results
# include("_batch_results/_batch_process_results.jl")
function p_tol_search(svm, pgrid, toldf, interpf; quietly=false)
    for i in range(1, stop=size(toldf, 1))
        # Find V satisfying the Slope condition:
        slope_cond = (abs.(interpf[:eq_deriv](pgrid)) .<= toldf[i, :eq_deriv])

        # Find V satisfying the Debt at Par Condition:
        aggP = [get_agg_p(svm, p=p) for p in pgrid]
        debt_at_par_cond = (abs.(interpf[:debt](pgrid) .- aggP) .<
                            toldf[i, :debt_diff])

        # Check if condition is satisfied:
        # Find Intersection of debt and equity conditions
        p_filtered = pgrid[.&(debt_at_par_cond, slope_cond)]

        if !isempty(p_filtered)
            if !quietly
                println("P Filter Conditions Satisfied! Exiting...")
            end
            return p_filtered
        end
    end
    return []
end


function p_interp_fun(svm, x::DataFrame, toldf::DataFrame; N::Integer=10^5, quietly::Bool=false)
    c = unique(x[:, :c])[1]
    pgrid = range(minimum(x[:, :p]), stop=maximum(x[:, :p]), length=N)

    interpf = Dict()
    for col in [:eq_deriv, :vb, :eq_min_val, :debt, :equity]
        interpf[col] = Dierckx.Spline1D(x[:, :p], x[:, col]; k=3, bc="extrapolate")
    end

    # Filter by (i) Debt Principal Difference,
    #           (ii) Equity Derivative
    p_filtered = p_tol_search(svm, pgrid, toldf,
                              interpf,
                              quietly=quietly)

    # Take the last occurrence of the minimum
    # that is, the largest VB value yielding the smallest derivative.
    p_filter_success = false
    if !isempty(p_filtered)
        p_filter_success = true
        if !quietly
            println(string("c: ", c, " -> P Filter success!"))
        end
        # inv_p_filtered = reverse(p_filtered)

        # Back-out solution -> Equity Derivative
        abs_debt_diffs = abs.(interpf[:debt](p_filtered) .-
                              [get_agg_p(svm, p=p) for p in p_filtered])
        optp = p_filtered[argmin(abs_debt_diffs)]
    else
        if !quietly
            println(string("c: ", c, " -> P Filter failed..."))
        end

        # Back-out solution -> Debt-At-Par + Equity Derivative:
        aggP = [get_agg_p(svm, p=p) for p in pgrid]
        debt_at_par_cond = reverse(abs.(interpf[:debt](pgrid) .- aggP))
        optp = reverse(pgrid)[argmin(debt_at_par_cond)]
        # eq_deriv_cond = abs.(interpf[:eq_deriv](inv_p_filtered))
        # optp = inv_p_filtered[argmin(.75 * debt_at_par_cond .+ .25 * eq_deriv_cond)]
    end

    # Back-out solutions:
    opt_debt = interpf[:debt](optp)
    aggP = get_agg_p(svm, p=optp)
    opt_eq = interpf[:equity](optp)
    opt_firm_val = opt_debt + opt_eq
    # return DataFrame(c = c,
    #                  p = optp,
    #                  opt_vb = interpf[:vb](optp),
    #                  cvml_vb = get_cvm_vb(svm, svm.pm.sigmal;
    #                                       mu_b=svm.mu_b, c=c, p=optp),
    #                  cvmh_vb = get_cvm_vb(svm, svm.pm.sigmah;
    #                                       mu_b=svm.mu_b, c=c, p=optp),
    #                  debt_diff = opt_debt - aggP,
    #                  debt_per_diff = (opt_debt - aggP) / aggP,
    #                  eq_deriv = interpf[:eq_deriv](optp),
    #                  eq_min_val = interpf[:eq_min_val](optp),
    #                  debt = opt_debt,
    #                  equity = opt_eq,
    #                  firm_value = opt_firm_val,
    #                  leverage = get_leverage(opt_debt, opt_eq),
    #                  MBR = get_mbr(svm.pm.V0, opt_debt, opt_eq),
    #                  p_filter_success = p_filter_success)
    return DataFrame(:c => c,
                     :p => optp,
                     :opt_vb => interpf[:vb](optp),
                     :cvml_vb => get_cvm_vb(svm, svm.pm.sigmal;
                                          mu_b=svm.mu_b, c=c, p=optp),
                     :cvmh_vb => get_cvm_vb(svm, svm.pm.sigmah;
                                          mu_b=svm.mu_b, c=c, p=optp),
                     :debt_diff => opt_debt - aggP,
                     :debt_per_diff => (opt_debt - aggP) / aggP,
                     :eq_deriv => interpf[:eq_deriv](optp),
                     :eq_min_val => interpf[:eq_min_val](optp),
                     :debt => opt_debt,
                     :equity => opt_eq,
                     :firm_value => opt_firm_val,
                     :leverage => get_leverage(opt_debt, opt_eq),
                     :MBR => get_mbr(svm.pm.V0, opt_debt, opt_eq),
                     :p_filter_success => p_filter_success)
    # leverage = (opt_debt / opt_firm_val) * 100,
    # ROE = (opt_eq / (svm.pm.V0 - opt_debt) - 1) * 100,

end


function process_combination_results(bt, svm;
                                     toldf::DataFrame=toldf,
                                     use_all_eqdf::Bool=true,
                                     drop_fail::Bool=false,
                                     save_df::Bool=true,
                                     dfname="soldf" )


    # Load Equity Finite Differences Files (eqdf_final)
    eqfds_final = [x for x in readdir(bt.mi.comb_res_path) if
                   .&(occursin("eq_fd", x), !occursin("all", x))]

    LL = @time fetch(@spawn [load_eq_results(bt, svm, dfn;
                                             use_all_eqdf=use_all_eqdf)
                             for dfn in eqfds_final])
    LL = [x for x in LL if !isempty(x)]
    eqdf_final = sort(vcat(LL...), [:c, :p])

    # Drop duplicates
    unique!(eqdf_final, [:c, :p])

    # For each coupon value, interpolate and extract results
    LL = @time fetch(@spawn [p_interp_fun(svm,
                             eqdf_final[abs.(eqdf_final[:, :c] .- c).<1e-4, :],
                             toldf) for c in unique(eqdf_final[:, :c])])
    soldf = sort(vcat(LL...), [:c])

    # Add Columns with Parameter Values
    cols = [x for x in vcat(bt.dfc.main_params, bt.dfc.fixed_params, [:mu_b, :m])
            if x !=:delta]
    for col in cols
        soldf[!, col] .= bt.mi._svm_dict[col]
    end
    soldf[!, :delta] .= soldf[:, :gross_delta] .- soldf[:, :iota]

    # Reoder columns
    cols1 = vcat(bt.dfc.main_params, [x for x in bt.dfc.k_struct_params if x !=:vb])
    cols2 = [x for x in Symbol.(names(soldf)) if .&(!(x in cols1), !(x in bt.dfc.fixed_params))]
    soldf = unique!(soldf[:, vcat(cols1, cols2, bt.dfc.fixed_params)])

    # Save DataFrame
    if save_df
        CSV.write(string(bt.mi.comb_res_path, "/", dfname, ".csv"), soldf)
    end

    return soldf
end


# ** Batch Filter
# include("_batch_results/_batch_filter.jl")
# Function below was adapted to run in Julia 1.1.0. It requires the packages
# 1. LinearAlgebra
# 2. DSP

#Polynomial smoothing with the Savitsky Golay filters
#
# Sources
# ---------
# Theory: http://www.ece.rutgers.edu/~orfanidi/intro2sp/orfanidis-i2sp.pdf
# Python Example: http://wiki.scipy.org/Cookbook/SavitzkyGolay
function savitsky_golay(x::Vector, windowSize::Integer, polyOrder::Integer; deriv::Integer=0)
    #Some error checking
    @assert isodd(windowSize) "Window size must be an odd integer."
    @assert polyOrder < windowSize "Polynomial order must me less than window size."

    halfWindow = Int((windowSize-1)/2)

    #Setup the S matrix of basis vectors.
    S = zeros(windowSize, polyOrder+1)
    for ct = 0:polyOrder
	#S[:,ct+1] = [-halfWindow:halfWindow].^(ct)
        S[:,ct+1] = range(-halfWindow, stop=halfWindow).^(ct)
    end

    #Compute the filter coefficients for all orders
    #From the scipy code it seems pinv(S) and taking rows should be enough
    G = S*LinearAlgebra.pinv(S'*S)

    #Slice out the derivative order we want
    filterCoeffs = G[:,deriv+1] * factorial(deriv)

    #Pad the signal with the endpoints and convolve with filter
    # paddedX = [x[1]*ones(halfWindow), x, x[end]*ones(halfWindow)]
    paddedX = vcat(x[1]*ones(halfWindow), x, x[end]*ones(halfWindow))
    # y = conv(filterCoeffs[end:-1:1], paddedX)
    y = DSP.conv(reverse(filterCoeffs), paddedX)

    #Return the valid midsection
    return y[2*halfWindow+1:end-2*halfWindow]

end


function findlocalmaxima(signal::Vector)
    inds = Int[]
    if length(signal)>1
        if signal[1]>signal[2]
            push!(inds,1)
        end

        for i=2:length(signal)-1
            if signal[i-1]<signal[i]>signal[i+1]
                push!(inds,i)
            end
        end

        if signal[end]>signal[end-1]
            push!(inds,length(signal))
        end
    end

    return inds
end


function filter_k_struct(df; interp_polyOrder::Integer=3,
                             filter_windowSize::Integer=5 * 10^2 + 1,
                         filter_polyOrder::Integer=3,
                         interp_vars=[:p, :opt_vb, :cvml_vb, :cvmh_vb, :debt, :equity])

    # Interpolate and Filter p, VB, Debt, Equity and Firm Value
    sgF = Dict()
    sgF[:cgrid] = range(minimum(df[:, :c]), stop=maximum(df[:, :c]), length=10^4)
    for x in interp_vars
        tmp = Dierckx.Spline1D(df[:, :c], df[:, x], k=interp_polyOrder, bc="extrapolate")
        sgF[x] = savitsky_golay(tmp(sgF[:cgrid]),  filter_windowSize, filter_polyOrder)
    end
    # sgF[:firm_value] = sgF[:debt] .+ sgF[:equity]

    return sgF
end


# ** Batch Optimal Capital Structure
# include("_batch_results/_batch_opt_k_struct.jl")
function opt_k_struct_cols(bt; de_cols::Array{Symbol,1}=[:debt_diff, :eq_deriv, :eq_min_val])
    share_cols = vcat([:sg_debt, :debt, :sg_equity],
                      [x for x in bt.dfc.share_values if x != :debt])
    debt_diff_cols = [x for x in bt.dfc.debt_vars if !(x in de_cols)]
    eqdf_cols = [x for x in bt.dfc.equity_vars if !(x in de_cols)]
    return vcat(:obj_fun, bt.dfc.main_params, bt.dfc.k_struct_params,
                [:cvml_vb], [:cvmh_vb],
                de_cols, share_cols, bt.dfc.fixed_params,
                debt_diff_cols, eqdf_cols)
end


function optimal_capital_struct(bt, svm;
                                firm_obj_fun::Symbol=:firm_value,
                                df::DataFrame=DataFrame(),
                                dfname::String=soldf_name,
                                interp_polyOrder::Integer=3,
                                filter_windowSize::Integer=5 * 10^2 + 1,
                                filter_polyOrder::Integer=3)

    # Load Solutions
    if isempty(df)
        df = CSV.read(string(bt.mi.comb_res_path, "/", dfname, ".csv"), DataFrame;
                      types=vcat(fill(Float64, 23), [Bool], fill(Float64, 5)))
    end

    # Discard failed results
    df = df[df[:, :p_filter_success].==true, :]

    # Interpolate and Filter Optimal Capital Structure
    sgF = filter_k_struct(df; interp_polyOrder=interp_polyOrder,
                          filter_windowSize=filter_windowSize,
                          filter_polyOrder=filter_polyOrder)

    # Compute Optimal Capital Structure
    sgF[:firm_value] = sgF[:debt] .+ sgF[:equity]
    sgF[:MBR] = get_mbr(get_param(svm, :V0), sgF[:debt], sgF[:equity])

    cpos = minimum(findlocalmaxima(sgF[firm_obj_fun]))

    # Run Equity Finite Differences Method
    eqDF = eq_fd(svm;
                 vbl=sgF[:opt_vb][cpos],
                 mu_b=svm.mu_b,
                 c=sgF[:cgrid][cpos],
                 p=sgF[:p][cpos])

    # Add Parameter Values
    eqDF[!, :cvml_vb] .= sgF[:cvml_vb][cpos]
    eqDF[!, :cvmh_vb] .= sgF[:cvmh_vb][cpos]
    eqDF[!, :sg_debt] .= sgF[:debt][cpos]
    eqDF[!, :sg_equity] .= sgF[:equity][cpos]

    # Add Debt Functions
    eqDF = debt_at_par_diffs(svm, eqDF, :p)

    # Add Column with firm_obj_fun
    eqDF[!, :obj_fun] .= String(firm_obj_fun)

    # Reoder columns
    cols = opt_k_struct_cols(bt)
    cols = [x for x in cols if (x in Symbol.(names(eqDF)))]

    return eqDF[:, cols]
end

function save_combination_opt_k_struct(bt, eqDF::DataFrame;
                                       idcols::Array{Symbol, 1}=vcat([:m, :obj_fun],
                                                                     bt.dfc.main_params,
                                                                     bt.dfc.fixed_params),
                                       dfname::String=optdf_name,
                                       coltypes::Array{DataType,1}=opt_k_struct_df_coltypes)

    if (string(dfname, ".csv") in readdir(bt.mi.comb_res_path))
        try
            seqDF = CSV.read(string(bt.mi.comb_res_path, "/", dfname, ".csv"), DataFrame;
                             types=coltypes)

            cond = .&(vcat([seqDF[:, :obj_fun] .== eqDF[:, :obj_fun]],
                           [(abs.(seqDF[:, col] .- eqDF[:, col]) .< 1e-6) for col in idcols if col != :obj_fun])...)
            if sum(cond) == 0
                println("No match found. Appending row...")
                seqDF = vcat(seqDF, eqDF)
            elseif sum(cond) == 1
                println("Match found! Replacing row...")
                seqDF[cond, :] .= eqDF
            else
                println("Multiple matches found. Refine ID columns. Returning...")
                return
            end
        catch
            println("Unable to load DataFrame.")
            seqDF = eqDF
        end
    else
        seqDF = eqDF
    end

    # Save DataFrame
    CSV.write(string(bt.mi.comb_res_path, "/", dfname, ".csv"), seqDF)
end


function compile_optimal_cap_struct(bt, svm;
                                    firm_obj_fun::Symbol=:firm_value,
                                    toldf::DataFrame=toldf,
                                    use_all_eqdf::Bool=true,
                                    drop_fail::Bool=false,
                                    save_soldf::Bool=true,
                                    soldf_name::String=soldf_name,
                                    interp_polyOrder::Integer=3,
                                    filter_windowSize::Integer=5 * 10^2 + 1,
                                    filter_polyOrder::Integer=3,
                                    replace_optdf::Bool=false,
                                    save_optdf::Bool=true,
                                    optdf_name::String=optdf_name,
                                    min_yield=1e-3,
                                    max_yield=5.,
                                    N=10^5,
                                    ftype::String="bf")
                  

    # Collect and Process Results
    soldf = process_combination_results(bt, svm;
                                        toldf=toldf,
                                        use_all_eqdf=use_all_eqdf,
                                        drop_fail=drop_fail,
                                        save_df=save_soldf,
                                        dfname=soldf_name)


    # Find Optimal Capital Structure
    eqDF = optimal_capital_struct(bt, svm;
                                  firm_obj_fun=firm_obj_fun,
                                  df=soldf,
                                  dfname=soldf_name,
                                  interp_polyOrder=interp_polyOrder,
                                  filter_windowSize=filter_windowSize,
                                  filter_polyOrder=filter_polyOrder)

    # Add combination IDs
    combdf = get_batch_comb_num(bt)
    eqDF[!, :comb_num] .= combdf[1, :comb_num]
    eqDF[!, :m_comb_num] .= combdf[1, :comb_num]
    id_cols = [:comb_num, :m, :m_comb_num]


    cols_order = vcat(id_cols, [x for x in Symbol.(names(eqDF)) if !(x in id_cols)])
    eqDF = eqDF[:, cols_order]


    # Compute Yield and Yield Spread
    ## Set capital structure
    KS = Dict{Symbol, Float64}()
    for x in [:mu_b, :m, :c, :p]
      KS[x] = eqDF[1, x]
    end
    KS[:vbl] = eqDF[1, :vb]

    eqDF[!, :yield] .= get_bond_yield(svm; KS=KS, min_yield=min_yield,
                                     max_yield=max_yield, N=N, ftype=ftype)
    eqDF[!, :yield_spd] .= get_bond_spread(svm; KS=KS, min_yield=min_yield,
                                          max_yield=max_yield, N=N, ftype=ftype)


    if replace_optdf
        CSV.write(string(bt.mi.comb_res_path, "/", optdf_name, ".csv"), eqDF)
    elseif save_optdf
        save_combination_opt_k_struct(bt, eqDF; dfname=optdf_name)
    end

    return eqDF
end


function get_opt_results(bt;
                         comb_num::Integer=0,
                         firm_obj_fun::Symbol=:firm_value,
                         use_all_eqdf::Bool=true,
                         m::Float64=0.,
                         m_comb_num::Integer=0,
                         display_msgs::Bool=true,
                         toldf::DataFrame=Batch.toldf,
                         drop_fail::Bool=false,
                         save_soldf::Bool=true,
                         soldf_name::String=soldf_name,
                         interp_polyOrder::Integer=3,
                         filter_windowSize::Integer=5 * 10^2 + 1,
                         filter_polyOrder::Integer=3,
                         replace_optdf::Bool=false,
                         save_optdf::Bool=true,
                         optdf_name::String=optdf_name)

    println(string("COMBNUM: ", comb_num))
    bt, svm = get_bt_svm(;comb_num=comb_num,
                         m=m, m_comb_num=m_comb_num,
                         display_msgs=display_msgs)

    try
        return compile_optimal_cap_struct(bt, svm;
                                          firm_obj_fun=firm_obj_fun,
                                          toldf=toldf,
                                          use_all_eqdf=use_all_eqdf,
                                          drop_fail=drop_fail,
                                          save_soldf=save_soldf,
                                          soldf_name=soldf_name,
                                          interp_polyOrder=interp_polyOrder,
                                          filter_windowSize=filter_windowSize,
                                          filter_polyOrder=filter_polyOrder,
                                          replace_optdf=replace_optdf,
                                          save_optdf=save_optdf,
                                          optdf_name=optdf_name)


        # return eqDF[vcat([:comb_num], [x for x in names(eqDF) if x !=:comb_num])]
    catch
        cols = opt_k_struct_cols(bt)
        eqDict = Dict()
        for col in cols
            if string(col) in keys(bt.mi._svm_dict)
               eqDict[col] = bt.mi._svm_dict[:, string(col)]
            elseif col != :eq_negative
               eqDict[col] = NaN
            else
               eqDict[col] = true
            end
        end
        eqDF = DataFrame(eqDict)
        # eqDF = hcat(get_batch_comb_num(bt)[[:comb_num, :m_comb_num]], eqDF)
        # id_cols = [:comb_num, :m, :m_comb_num]
        # cols_order = vcat(id_cols, [x for x in names(eqDF) if !(x in id_cols)])
        # return eqDF[cols_order]
        # Add combination IDs
        eqDF[!, :comb_num] .= comb_num
        combdf = get_batch_comb_num(bt)
        optDF[!, :comb_num] .= combdf[1, :comb_num]
        optDF[!, :m_comb_num] .= combdf[1, :comb_num]
        id_cols = [:comb_num, :m, :m_comb_num]
        cols_order = vcat(id_cols, [x for x in names(eqDF) if !(x in id_cols)])
        return eqDF[:, cols_order]
    end
end


function compile_svm_opt_results(bt; m::Float64=NaN,
                                 firm_obj_fun::Symbol=:firm_value,
                                 use_all_eqdf::Bool=true,
                                 display_msgs::Bool=false,
                                 toldf::DataFrame=Batch.toldf,
                                 drop_fail::Bool=false,
                                 save_soldf::Bool=true,
                                 soldf_name::String=soldf_name,
                                 interp_polyOrder::Integer=3,
                                 filter_windowSize::Integer=5 * 10^2 + 1,
                                 filter_polyOrder::Integer=3,
                                 save_optdf::Bool=true,
                                 optdf_name::String=optdf_name,
                                 recompute_comb_opt_res::Bool=true,
                                 replace_optdf::Bool=false,
                                 save_results::Bool=true,
                                 opt_k_struct_df_name::String=opt_k_struct_df_name)

    if !isnan(m)
        comb_nums = bt.bp.df[bt.bp.df[:, :m] .== m, :comb_num]
    else
        comb_nums = bt.bp.df[:, :comb_num]
    end

    #    comb_nums = range(1, stop=size(hcat(bt._params_combs...)', 1))
    if recompute_comb_opt_res
        optDF_LL = @time fetch(@spawn [get_opt_results(bt;
                                                       comb_num=comb_num,
                                                       firm_obj_fun=firm_obj_fun,
                                                       use_all_eqdf=use_all_eqdf,
                                                       display_msgs=display_msgs,
                                                       toldf=toldf,
                                                       drop_fail=drop_fail,
                                                       save_soldf=save_soldf,
                                                       soldf_name=soldf_name,
                                                       interp_polyOrder=interp_polyOrder,
                                                       filter_windowSize=filter_windowSize,
                                                       filter_polyOrder=filter_polyOrder,
                                                       replace_optdf=replace_optdf,
                                                       save_optdf=save_optdf,
                                                       optdf_name=optdf_name)
                                       for comb_num in comb_nums])
    else
        optDF_LL = @time fetch(@spawn [load_svm_opt_results(bt;
                                                            comb_num=comb_num,
                                                            display_msgs=display_msgs,
                                                            toldf=toldf,
                                                            drop_fail=drop_fail,
                                                            save_soldf=save_soldf,
                                                            soldf_name=soldf_name,
                                                            interp_polyOrder=interp_polyOrder,
                                                            filter_windowSize=filter_windowSize,
                                                            filter_polyOrder=filter_polyOrder,
                                                            save_optdf=save_optdf,
                                                            optdf_name=optdf_name)
                                       for comb_num in comb_nums])
    end
    optDF = vcat(optDF_LL...)

    # Filter
    optDF = optDF[optDF[:, :obj_fun] .== String(firm_obj_fun), :]

    if save_results
        println("Saving compiled results...")
        extension = lowercase(string(firm_obj_fun))
        dfname = string(opt_k_struct_df_name, "_", extension)

        if isnan(m)
            bt = set_par_dict(bt; comb_num=1, display_msgs=display_msgs)
            bt = set_comb_res_paths(bt)
            fpath = bt.mi.batch_res_path
        else
            bt = set_par_dict(bt; m=m, m_comb_num=1, display_msgs=display_msgs)
            bt = set_comb_res_paths(bt)
            fpath = bt.mi.maturity_path
        end
        CSV.write(string(fpath, "/", dfname, ".csv"), optDF)
    end

    return optDF
end


# * Batch CVM Methods
# include("_batch_cvm_funs.jl")
function get_k_struct(cvm, mu_b::Float64, m::Float64, c::Float64, p::Float64)
    if isnan(mu_b)
        mu_b = cvm.mu_b
    end
    if isnan(m)
        m = cvm.m
    end
    if isnan(c)
        c = cvm.c
    end
    if isnan(p)
        p = cvm.p
    end
    return mu_b, m, c, p
end


function get_p_grid(cvm; mu_b::Float64=NaN, m::Float64=NaN,
                    c::Float64=NaN, step::Float64=.2, N1::Int64=10^2, N2::Int64=10^4)
    mu_b, m, c, _ = get_k_struct(cvm, mu_b, m, c, NaN)
    rdisc = get_rdisc(cvm)
    C = get_agg_c(cvm; mu_b=mu_b, c=c)

    # Create grid for P:
    # 1. Set Max p:
    if (m - 1/rdisc) * (1 - exp(rdisc * m)) > 0
        pMax = (C/ m) / rdisc
    else
        pMax = 3.0*cvm.pm.V0 / m
    end

    # 2. Search for minimum acceptable P value. This is the value of P for which vb = 0:
    p_grid = range(.5, stop=cvm.pm.V0, length=N1)
    vb_grid = fetch(@spawn [get_cvm_vb(cvm, cvm.pm.sigmal;
                                       mu_b=mu_b, c=c, p=p) for p in p_grid])
    vbff = Dierckx.Spline1D(p_grid, vb_grid, k=3, bc="extrapolate")
    p_grid_ref = range(p_grid[1], stop=p_grid[end], length=N2)

    # vb must be positive
    if sum(vbff(p_grid_ref) .< 0) > 0
        p_lower_bound = maximum(p_grid_ref[(vbff(p_grid_ref) .< 0)])
        p_grid_ref = p_grid_ref[p_grid_ref .> p_lower_bound]
    end
    p0 = p_grid_ref[argmin(abs.(vbff(p_grid_ref) .- 1e-4))]

    # 3. Return Grid (start at 1.1*p0, to avoid potential function errors)
    return range(1.001 * p0, stop = 1.1 * pMax, step = step)
end


function debt_price_diff(cvm, mu_b::Float64, c::Float64, p::Float64)
    vbl = get_cvm_vb(cvm, cvm.pm.sigmal; mu_b=mu_b, c=c, p=p)
    debt = get_cvm_debt_price(cvm, vbl, cvm.pm.sigmal; mu_b=mu_b, c=c, p=p)

    return debt - get_agg_p(cvm, mu_b=mu_b, p=p)
end


function get_p_debt_at_par(cvm; mu_b::Float64=NaN, m::Float64=NaN, c::Float64=NaN,
                           step::Float64=.2, N1::Int64=10^2, N2::Int64=10^4)
    mu_b, m, c, _ = get_k_struct(cvm, mu_b, m, c, NaN)

    # Find P such that Debt = P

    # Given C, we search for P such that P = D, at optimal vb.
    # The function below can have more than one root, so initial
    # guess matters!
    # The solution is to the guess the lowest p value possible: 0
    # To avoid this issue, I rely on grid search.

    # Set p grid:
    pGrid = get_p_grid(cvm; mu_b=mu_b, m=m, c=c, step=step, N1=N1, N2=N2)

    # Compute squared difference between Debt and P:
    pDebtDiff = fetch(@spawn [debt_price_diff(cvm, mu_b, c, p) for p in pGrid])

    diff = Dierckx.Spline1D(pGrid, pDebtDiff, k=3, bc="extrapolate")
    pvec = range(pGrid[1], stop=pGrid[end], length=N2)

    return  pvec[argmin(abs.(diff(pvec)))]
end


function get_cvm_c_debt_at_par(bt, cvm, c::Float64;
                               mu_b::Float64=NaN, m::Float64=NaN,
                               step::Float64=.2, N1::Int64=10^2, N2::Int64=10^4)
    # Capital Structure
    mu_b, m, c, _ = get_k_struct(cvm, mu_b, m, c, NaN)
    # P s.t. P = Debt
    pOpt = get_p_debt_at_par(cvm; mu_b=mu_b, m=m, c=c, step=step, N1=N1, N2=N2)

    # Store Results ##############################################

    # Parameters
    varlist = [var for var in vcat(bt.dfc.main_params, bt.dfc.fixed_params) if var != :delta]
    tmp = Dict()
    for var in varlist
        tmp[var] = bt.mi._svm_dict[var]
    end
    df = DataFrame(tmp)

    df[!, :delta] .= df[:, :gross_delta] - df[:, :iota]

    df[!, :mu_b] .= mu_b
    df[!, :m] .= m
    df[!, :c] .= c
    df[!, :p] .= pOpt

    # Default Barrier
    vbl = get_cvm_vb(cvm, cvm.pm.sigmal; mu_b=mu_b, c=c, p=pOpt)
    df[!, :vb] .= vbl

    # Debt
    df[!, :debt] .= get_cvm_debt_price(cvm, vbl, cvm.pm.sigmal; mu_b=mu_b, c=c, p=pOpt)
    df = debt_at_par_diffs(cvm, df, :p)

    # Equity
    df[!, :eq_vb] .= get_cvm_eq(cvm, vbl, cvm.pm.sigmal; mu_b=mu_b, c=c, p=pOpt)
    df[!, :eq_min_val] .= df[1, :eq_vb] #[1]
    df[!, :equity] .= get_cvm_eq(cvm, cvm.pm.V0, cvm.pm.sigmal; mu_b=mu_b, c=c, p=pOpt)

    # Firm Value, Leverage & MBR
    df = non_interp_values(cvm, df)

    # Rearrange Columns
    return df[:, bt.dfc.dfcols]
end


function get_cvm_debt_at_par(bt, cvm;
                             mu_b::Float64=NaN, m::Float64=NaN,
                             step::Float64=.2, N1::Int64=10^2, N2::Int64=10^4,
                             save_soldf::Bool=true,
                             soldf_name::String=soldf_name)

    resdfs = fetch(@spawn [get_cvm_c_debt_at_par(bt, cvm, c;
                                                 mu_b=mu_b, m=m,
                                                 step=step, N1=N1, N2=N2)
                           for c in bt.coupon_grid])

    soldf = vcat(resdfs...)

    if save_soldf
        println("Saving CVM results...")
        CSV.write(string(bt.mi.comb_res_path, "/", soldf_name, ".csv"), soldf)
    end

    return soldf
end


function optimal_cvm_capital_struct(bt, cvm;
                                    firm_obj_fun::Symbol=:firm_value,
                                    df::DataFrame=DataFrame(),
                                    dfname::String=Batch.soldf_name,
                                    interp_polyOrder::Int64=3,
                                    filter_windowSize::Int64=5 * 10^2 + 1,
                                    filter_polyOrder::Int64=3,
                                    replace_optdf::Bool=false,
                                    save_optdf::Bool=true,
                                    optdf_name::String=optdf_name)

    # Load Solutions
    if isempty(df)
        df = CSV.read(string(bt.mi.comb_res_path, "/", dfname, ".csv"), DataFrame;
                      types=fill(Float64, 28))
    end

    # Interpolate and Filter Optimal Capital Structure
    sgF = filter_k_struct(df; interp_polyOrder=interp_polyOrder,
                          filter_windowSize=filter_windowSize,
                          filter_polyOrder=filter_polyOrder,
                          interp_vars=[:p, :vb, :debt, :equity])

    # Compute Optimal Capital Structure
    sgF[:firm_value] = sgF[:debt] .+ sgF[:equity]
    sgF[:MBR] = get_mbr(get_param(cvm, :V0), sgF[:debt], sgF[:equity])


    # Compute Optimal Capital Structure
    cpos = minimum(findlocalmaxima(sgF[firm_obj_fun]))
    optDF = DataFrame()

    optDict = Dict()
    # Parameters ######################
    for var in vcat(bt.dfc.main_params, bt.dfc.fixed_params, [:mu_b, :m])
        optDict[var] = df[1, var]
    end
    # #################################


    #  Capital Struct #################
    optDict[:c] = sgF[:cgrid][cpos]
    optDict[:p] = sgF[:p][cpos]
    optDict[:vb] = sgF[:vb][cpos]
    # #################################

    # Debt ###########################
    optDict[:sg_debt] = sgF[:debt][cpos]

    optDict[:debt] = get_cvm_debt_price(cvm, optDict[:vb], #[1],
                                        optDict[:sigmal]; #[1];
                                        mu_b=optDict[:mu_b], #[1],
                                        c=optDict[:c], # [1],
                                        p=optDict[:p]) # [1])
    optDF = DataFrame(optDict)
    optDF = debt_at_par_diffs(cvm, optDF, :p)
    # #################################

    # Equity #########################
    optDF[!, :eq_vb] .= get_cvm_eq(cvm,
                                  optDF[1, :vb], #[1],
                                  optDF[1, :sigmal]; #[1];
                                  mu_b=optDF[1, :mu_b], #[1],
                                  c=optDF[1, :c], #[1],
                                  p=optDF[1, :p]) #[1])
    optDF[!, :eq_min_val] .= optDF[1, :eq_vb][1]
    optDF[!, :sg_equity] .= sgF[:equity][cpos]
    optDF[!, :equity] .= get_cvm_eq(cvm,
                                   optDF[1, :V0], #[1],
                                   optDF[1, :sigmal]; #[1];
                                   mu_b=optDF[1, :mu_b], #[1],
                                   c=optDF[1, :c], #[1],
                                   p=optDF[1, :p]) #[1])
    # #################################
    optDF = non_interp_values(cvm, optDF)

    # CVM NaN columns
    optDF[!, :cvml_vb] .= NaN
    optDF[!, :cvmh_vb] .= NaN
    optDF[!, :eq_deriv] .= NaN
    optDF[!, :eq_min_val] .= NaN
    optDF[!, :eq_deriv_min_val] .= NaN

    optDF[!, :eq_negative] .= false

    # Add Column with firm_obj_fun
    optDF[!, :obj_fun] .= String(firm_obj_fun)

    # Add combination IDs
    combdf = get_batch_comb_num(bt)
    optDF[!, :comb_num] .= combdf[1, :comb_num]
    optDF[!, :m_comb_num] .= combdf[1, :comb_num]
    id_cols = [:comb_num, :m, :m_comb_num]
    cols_order = vcat(id_cols, [x for x in opt_k_struct_cols(bt) if !(x in id_cols)])
    optDF = optDF[:, cols_order]


    if replace_optdf
        CSV.write(string(bt.mi.comb_res_path, "/", optdf_name, ".csv"), optDF)
    elseif save_optdf
        idcols = [x for x in vcat([:m, :obj_fun],
                                  bt.dfc.main_params,
                                  bt.dfc.fixed_params) if !(x in [:lambda, :sigmah])]
        save_combination_opt_k_struct(bt, optDF; dfname=optdf_name,
                                      idcols=idcols,
                                      coltypes=cvm_opt_k_struct_df_coltypes)
    end

    return optDF
end


function cvm_results_loader(comb_num::Int64;
                            optdf_name::String=optdf_name,
                            coltypes::Array{DataType,1}=cvm_opt_k_struct_df_coltypes,
                            display_msgs::Bool=false)

    bt = get_bt(; model="cvm", comb_num=comb_num, display_msgs=display_msgs)
    optDF = CSV.read(string(bt.mi.comb_res_path, "/", optdf_name, ".csv"), DataFrame; types=coltypes)

    return optDF
end


function get_cvm_opt_results(; comb_num::Integer=0,
                             firm_obj_fun::Symbol=:firm_value,
                             m::Float64=0.,
                             m_comb_num::Integer=0,
                             display_msgs::Bool=true,
                             dfname::String=soldf_name,
                             replace_optdf::Bool=false,
                             save_optdf::Bool=true,
                             optdf_name::String=optdf_name)

    bt, cvm = get_bt_cvm(;comb_num=comb_num,
                         m=m, m_comb_num=m_comb_num,
                         display_msgs=display_msgs)

    return optimal_cvm_capital_struct(bt, cvm;
                                      firm_obj_fun=firm_obj_fun,
                                      dfname=dfname,
                                      replace_optdf=replace_optdf,
                                      save_optdf=save_optdf,
                                      optdf_name=optdf_name)
end


function compile_cvm_opt_results(; m::Float64=NaN,
                                 firm_obj_fun::Symbol=:firm_value,
                                 display_msgs::Bool=false,
                                 dfname::String=soldf_name,
                                 replace_optdf::Bool=false,
                                 save_optdf::Bool=true,
                                 optdf_name::String=optdf_name,
                                 save_results::Bool=true,
                                 opt_k_struct_df_name::String=opt_k_struct_df_name)

    bt = get_bt(; model="cvm", comb_num=1, display_msgs=display_msgs)

    if !isnan(m)
        comb_nums = bt.bp.df[abs.(bt.bp.df[:, :m] .- m) .< 1e-6, :comb_num]
    else
        comb_nums = bt.bp.df[:, :comb_num]
    end

    optDFs = fetch(@spawn [get_cvm_opt_results(; comb_num=comb,
                                               firm_obj_fun=firm_obj_fun,
                                               display_msgs=display_msgs,
                                               dfname=dfname,
                                               replace_optdf=replace_optdf,
                                               save_optdf=save_optdf,
                                               optdf_name=optdf_name)
                           for comb in comb_nums])

    cvmOptDF = sort!(vcat(optDFs...), [:m, :m_comb_num])

    if save_results
        println("Saving compiled results...")
        extension = lowercase(string(firm_obj_fun))
        dfname = string(opt_k_struct_df_name, "_", extension)

        if isnan(m)
            bt = set_par_dict(bt; comb_num=1, display_msgs=display_msgs)
            bt = set_comb_res_paths(bt)
            fpath = bt.mi.batch_res_path
        else
            bt = set_par_dict(bt; m=m, m_comb_num=1, display_msgs=display_msgs)
            bt = set_comb_res_paths(bt)
            fpath = bt.mi.maturity_path
        end
        CSV.write(string(fpath, "/", dfname, ".csv"), cvmOptDF)
    end

    return cvmOptDF
end


function load_cvm_opt_results_df(; m::Float64=NaN,
                                 optdf_name::String=opt_k_struct_df_name,
                                 firm_obj_fun::Symbol=:firm_value,
                                 coltypes::Array{DataType,1}=cvm_opt_k_struct_df_coltypes,
                                 display_msgs::Bool=false,
                                 opt_k_struct_df_name::String=opt_k_struct_df_name)

    bt = get_bt(; model="cvm", comb_num=1, display_msgs=display_msgs)

    extension = lowercase(string(firm_obj_fun))
    dfname = string(opt_k_struct_df_name, "_", extension)
    if isnan(m)
        bt = set_par_dict(bt; comb_num=1, display_msgs=display_msgs)
        bt = set_comb_res_paths(bt)
        fpath = bt.mi.batch_res_path
    else
        bt = set_par_dict(bt; m=m, m_comb_num=1, display_msgs=display_msgs)
        bt = set_comb_res_paths(bt)
        fpath = bt.mi.maturity_path
    end

    try
        println("Loading optimal results dataframe...")
        return CSV.read(string(fpath, "/", dfname, ".csv"), DataFrame; types=cvm_opt_k_struct_df_coltypes)
    catch
        println("Unable to load optimal results dataframe. Recomputing...")
        return compile_cvm_opt_results(; m=m,
                                       firm_obj_fun=firm_obj_fun,
                                       opt_k_struct_df_name=opt_k_struct_df_name)
    end
end


# * END MODULE
end
