module_path = "/home/artur/BondPricing/Julia/modules/"
modnames = ["ModelObj", "AnalyticFunctions", "BondPrInterp", "EqFinDiff"]
for modl in modnames
    if !(joinpath(module_path, modl) in LOAD_PATH)
        push!(LOAD_PATH, joinpath(module_path, modl))
    end
end


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
using AnalyticFunctions: get_cvm_vb, get_agg_c, get_agg_p, get_rdisc
using BondPrInterp: get_cvm_debt_price, get_svm_debt_price, bpr_surfs, bpr_interp
using EqFinDiff: get_cvm_eq, eq_fd


include("_batch_inputs.jl")
include("_batch_structs.jl")

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

# Batch Jobs - Set & Get Methods
include("_batch_set_get_methods/_batch_set_get_params.jl")
include("_batch_set_get_methods/_batch_set_get_dirs.jl")
include("_batch_set_get_methods/_batch_get_objs.jl")

# Batch Jobs - File Methods
include("_batch_file_methods/_batch_obj_files.jl")
include("_batch_file_methods/_batch_equity_files.jl")
include("_batch_file_methods/_batch_results_files.jl")

# Batch Jobs - Debt at Par and Equity Finite Differences Method
include("_batch_core/_batch_pvb_grids.jl")
include("_batch_core/_batch_debt_at_par_funs.jl")
include("_batch_core/_batch_eq_fd.jl")

# Processing Results
include("_batch_results/_batch_diagnosis.jl")
include("_batch_results/_batch_process_results.jl")
include("_batch_results/_batch_filter.jl")
include("_batch_results/_batch_opt_k_struct.jl")

# CVM Methods
include("_batch_cvm_funs.jl")

end

