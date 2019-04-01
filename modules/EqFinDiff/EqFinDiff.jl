module_path = "/home/artur/BondPricing/Julia/modules/"
modnames = ["ModelObj", "AnalyticFunctions", "BondPrInterp"]
for modl in modnames
    if !(joinpath(module_path, modl) in LOAD_PATH)
        push!(LOAD_PATH, joinpath(module_path, modl))
    end
end


module EqFinDiff

using Distributed # @spawn macro
using Interpolations
using LinearAlgebra
using DataFrames
using Dierckx

# using PyCall
# @pyimport scipy.interpolate as sp_interpolate

# I am using functions from the ANalytic Functions module
# rfbond_price

# User-defined package:
using ModelObj: get_obj_model

using AnalyticFunctions: rdisc_pvs, cvm_eq,
                         get_rgrow, get_rdisc,
                         get_cvm_vb, get_param,
                         get_k_struct,
                         get_agg_c, get_agg_p,
                         on_default_payoff, rfbond_price

using BondPrInterp: grid_creator,
                    get_pv_rfdebt,
                    get_cvm_bond_price,
                    get_cvm_debt_price,
                    get_svm_bond_price,
		    get_svm_debt_price

include("eq_get_methods.jl")
include("eq_fin_diff_funs.jl")
# include("tmp_file.jl")

end


# fun_rel_path  = "../functions/"
# modname = "EqFinDiff"
# include(joinpath(fun_rel_path, modname, "eq_get_methods.jl")) 
# include(joinpath(fun_rel_path, modname, "eq_fin_diff_funs.jl")) 
# include(joinpath(fun_rel_path, modname, "tmp_file.jl")) 

