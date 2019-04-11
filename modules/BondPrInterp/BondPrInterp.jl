module_path = "/home/artur/BondPricing/Julia/modules/"
# modlpath = joinpath(module_path, "AnalyticFunctions")
# if !(modlpath in LOAD_PATH)
#     push!(LOAD_PATH, modlpath)
# end

modnames = ["ModelObj", "AnalyticFunctions"]
for modl in modnames
    if !(joinpath(module_path, modl) in LOAD_PATH)
        push!(LOAD_PATH, joinpath(module_path, modl))
    end
end



# Interpolate Primitive SVM Bond Pricing Functions
module BondPrInterp

using Distributed  # @spawn macro
using Distributions
using Interpolations
# using PyCall
# @pyimport scipy.interpolate as sp_interpolate
# @pyimport numpy as np
# using LinearAlgebra
using Dierckx

# User-defined packages:
using ModelObj: grid_creator, set_bpr_grids
using AnalyticFunctions: rdisc, rdisc_pvs, rfbond_price,
                         cvm_bond_price, rf_debt,
                         get_rdisc, get_rgrow, get_k_struct,
                         get_cvm_vb, get_param,
                         psi_v_td, dv_psi_v_td,
                         cvm_F, cvm_G, cvm_vb, get_cvm_vb,
                         on_default_payoff, no_vol_shock_cf_pv

include("svm_primitive_bpr_funs.jl")
include("svm_interp_bpr_funs.jl")
include("svm_bpr_fixed_ttm_funs.jl")
include("bondpr_funs.jl") 
include("debtpr_funs.jl") 

end